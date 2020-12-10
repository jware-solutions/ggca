// TODO: this file should be replaced when https://github.com/appaquet/extsort-rs/issues/5 is merged

use itertools::{Itertools, KMerge};
use rayon::prelude::*;
use tempdir::TempDir;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Result, Write};
use std::marker::PhantomData;
use std::path::{Path, PathBuf};

pub trait Sortable: Eq + Ord + Sized + Send + Sync {
    fn encode<W: Write>(&self, writer: &mut W);
    fn decode<R: Read>(reader: &mut R) -> Option<Self>;
}

pub struct ExternalSorter {
    max_mem: usize, // bytes
    tempdir: TempDir,
    file_count: usize,
}

impl ExternalSorter {
    pub fn new(max_mem: usize) -> ExternalSorter {
        assert!(max_mem > 0);
        ExternalSorter {
            max_mem, 
            tempdir: TempDir::new("sort_fasta").unwrap(),
            file_count: 0,
        }
    }

    fn file_name(&self, n: usize) -> PathBuf {
        self.tempdir.path().join(format!("{}", n))
    }

    fn write_chunk<T: Sortable>(&mut self, buf: &mut Vec<T>) -> Result<()> {
        buf.par_sort_unstable();
        let file = File::create(self.file_name(self.file_count))?;
        let mut out = BufWriter::new(file);
        for item in buf.iter() {
            item.encode(&mut out);
        }
        out.flush()?;
        self.file_count += 1;
        Ok(())
    }

    /// Sort a given iterator, returning a new iterator with items
    pub fn sort<T, I>(
        &mut self,
        iterator: I,
    ) -> Result<KMerge<<ChunkReader<T> as IntoIterator>::IntoIter>>
    where
        T: Sortable,
        I: Iterator<Item = T>,
    {
        let max_size = self.max_mem / std::mem::size_of::<T>() as usize;
        let mut buf = Vec::<T>::with_capacity(max_size);
        for item in iterator {
            if buf.len() == max_size {
                self.write_chunk(&mut buf)?;
                buf.clear();
            }
            buf.push(item);
        }
        if !buf.is_empty() {
            self.write_chunk(&mut buf)?;
        }
        std::mem::drop(buf);
        let mut readers = vec![];
        for n in 0..self.file_count {
            readers.push(ChunkReader::new(self.file_name(n))?);
        }
        Ok(readers.into_iter().kmerge())
    }
}

pub struct ChunkReader<T: Sortable> {
    reader: BufReader<std::fs::File>,
    phantom: PhantomData<T>,
}

impl<T: Sortable> ChunkReader<T> {
    pub fn new<P: AsRef<Path>>(fname: P) -> Result<ChunkReader<T>> {
        Ok(ChunkReader {
            reader: BufReader::new(File::open(fname)?),
            phantom: PhantomData,
        })
    }
}

impl<T: Sortable> Iterator for ChunkReader<T> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        T::decode(&mut self.reader)
    }
}
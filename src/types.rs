use crate::correlation::CorResult;

pub type VecOfResults = Vec<CorResult>;
// A tuple with Gene/GEM, Cpg Site ID (optional), and a vec of values
pub type TupleExpressionValues = (String, Option<String>, Vec<f64>);
pub type LazyMatrixInner = Box<dyn Iterator<Item = TupleExpressionValues>>;
pub type CollectedMatrix = Vec<TupleExpressionValues>;

import ggca


def main():
    mrna_file_path = "mrna.csv"
    gem_file_path = "mirna.csv"

    try:
        (result_combinations, _total_combinations_count, evaluated_combinations) = ggca.correlate(
            mrna_file_path,
            gem_file_path,
            correlation_method=ggca.CorrelationMethod.Pearson,
            correlation_threshold=0.5,
            sort_buf_size=2_000_000,
            adjustment_method=ggca.AdjustmentMethod.BenjaminiHochberg,
            is_all_vs_all=True,
            gem_contains_cpg=False,
            collect_gem_dataset=None,
            keep_top_n=2  # Keeps only top 2 elements
        )

        print(f'Number of resulting combinations: {len(result_combinations)} of {evaluated_combinations} evaluated '
              f'combinations')
        for combination in result_combinations:
            print(combination.gene, combination.gem, combination.correlation, combination.p_value,
                  combination.adjusted_p_value)
    except ggca.GGCADiffSamplesLength as ex:
        print('Raised GGCADiffSamplesLength:', ex)
    except ggca.GGCADiffSamples as ex:
        print('Raised GGCADiffSamples:', ex)
    except ggca.InvalidCorrelationMethod as ex:
        print('Raised InvalidCorrelationMethod:', ex)
    except ggca.InvalidAdjustmentMethod as ex:
        print('Raised InvalidAdjustmentMethod:', ex)
    except ggca.GGCAError as ex:
        print('Raised GGCAError:', ex)


if __name__ == '__main__':
    main()

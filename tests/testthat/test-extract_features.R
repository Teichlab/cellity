## Testing functions for extracting features ##

context("tests on feature extraction")

test_that("the extract_features() function works on example data", {
    data(sample_counts)
    data(sample_stats)
    sample_counts_nm <- normalise_by_factor(sample_counts, colSums(sample_counts))
    sample_features <- extract_features(sample_counts_nm, sample_stats)
    
    expect_that(sample_features, is_a("list"))
})



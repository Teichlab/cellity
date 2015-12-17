## Testing functions for extracting features ##

context("tests on feature extraction")

### Look at examples below
# 
# test_that("tests for presence of cellData or counts", {
#     expect_that(newSCESet(), 
#                 throws_error("Require at least one of exprsData, tpmData, fpkmData or countData argument."))
# })
# 
# test_that("example datasets work", {
#     data("sc_example_counts")
#     data("sc_example_cell_info")
#     pd <- new("AnnotatedDataFrame", data=sc_example_cell_info)
#     example_sceset <- newSCESet(countData=sc_example_counts, phenoData=pd)
#     example_sceset
#     
#     expect_that(example_sceset, is_a("SCESet"))
# })

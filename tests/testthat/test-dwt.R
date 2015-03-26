context("DWT - Unit Tests")


set.seed(999)
x = rnorm(2^8)


dwt_cpp(x, "haar", 4, boundary="periodic")


expect_equal(dwt_cpp(x, "haar", 4, boundary="periodic"))

expect_equal(qmf(g, inverse=TRUE), h )

expect_equal(qmf(g, inverse=FALSE), as.matrix(rev(h)) )


expect_equal(haar_filter(), wf.haar, check.attributes = FALSE)

expect_equal(select_filter(filter_name = "haar"), wf.haar, check.attributes = FALSE)

context("Wave Filters - Unit Tests")

g = as.matrix(rep(1/sqrt(2),2))

h = as.matrix(c(1/sqrt(2), -1/sqrt(2)))

wf.haar = list( as.matrix(2.0), #L
                h,              #h
                g)              #g

expect_equal(qmf(g), h )

expect_equal(qmf(g, inverse=TRUE), h )

expect_equal(qmf(g, inverse=FALSE), as.matrix(rev(h)) )


expect_equal(haar_filter(), wf.haar, check.attributes = FALSE)

expect_equal(select_filter(filter_name = "haar"), wf.haar, check.attributes = FALSE)

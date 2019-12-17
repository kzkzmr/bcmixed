library(testthat)
library(bcmixed)

context("Test bcmmrm")
data(aidscd4)
se1 <- bcmmrm(cd4, treatment, aidscd4, weekc, id,c("cd4.bl", "age", "sex"),
              c(0, 0, 1), "AR(1)")$meddif.rob.adj[[4]][2, 4]
se2 <- bcmmrm(cd4, treatment,
              aidscd4[aidscd4$weekc == 32, ])$meddif.rob.adj[[1]][2, 4]
test_that("Test bcmmrm not cleared", {
  expect_that(round(se1, 3), equals(1.292))
  expect_that(round(se2, 3), equals(1.718))
})

source("00_requirements.R")
source("01_GenerateData.R")
source("02_TTests.R")
source("03_Evaluate.R")

set.seed(29305)
dat <- generate_data()

# give both methods the readings and let farmtest estimate the factor count
alpha <- 0.05
farm <- FarmTest::farm.test(dat$X, KX = -1, alpha = alpha, p.method = "normal")
ordinary <- fit_t_tests(dat$X, alpha)
scores <- rbind(
  score_fit("t tests + BH", ordinary, dat$mu),
  score_fit("FarmTest", farm, dat$mu)
)

dir.create("results", showWarnings = FALSE)
write.csv(scores, "results/illustration_scores.csv", row.names = FALSE)
print(scores)

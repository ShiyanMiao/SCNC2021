
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("usethis")

library(usethis)
use_git_config(user.name = "Shiyan Miao", user.email = "miaoshiyan.annabel@gmail.com")
create_github_token()
gitcreds::gitcreds_set()

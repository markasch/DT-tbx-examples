# MC Integration (R version)

N <- 1000000
count <- 0

for (i in 1:N) {
  x <- runif(1)
  y <- runif(1)
  if (x^2 + y^2 < 1){
    count <- count + 1
  }
}

4*count/N

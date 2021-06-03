
## code to prepare `DATASET` dataset goes here
test_data <- array(0, dim = c(50, 50, 50))
for (x in 1:50){
    for (y in 1:50) {
        for (z in 1:50) {
             if (sqrt((x-10.1)^2+(y-35.1)^2 + (z - 35.1)^2) <= 7.5) {
                test_data[x, y, z] <- (1000 / sqrt((x-10.1)^2+(y-35.1)^2 + (z - 35.1)^2))
             }
        }
    }
}
tomoseq_test_x <- apply(test_data, 1, sum)
tomoseq_test_y <- apply(test_data, 2, sum)
tomoseq_test_z <- apply(test_data, 3, sum)
tomoseq_test_x <- t(tomoseq_test_x) %>% as.data.frame()
tomoseq_test_y <- t(tomoseq_test_y) %>% as.data.frame()
tomoseq_test_z <- t(tomoseq_test_z) %>% as.data.frame()

add_to_data <- function (data) {
    add_data <- runif(50) * 1000
    return(rbind(data, add_data))
}

for (i in 1:49) {
    tomoseq_test_x <- add_to_data(tomoseq_test_x)
    tomoseq_test_y <- add_to_data(tomoseq_test_y)
    tomoseq_test_z <- add_to_data(tomoseq_test_z)
}

row_name <- c()
for (i in 1:50) {
    row_name <- append(row_name, paste("gene", i, sep=""))
}
col_name <- as.character(c("gene_ID", seq(1,50)))
tomoseq_test_x <- cbind(row_name, tomoseq_test_x)
tomoseq_test_y <- cbind(row_name, tomoseq_test_y)
tomoseq_test_z <- cbind(row_name, tomoseq_test_z)

colnames(tomoseq_test_x) <- col_name
colnames(tomoseq_test_y) <- col_name
colnames(tomoseq_test_z) <- col_name
usethis::use_data(tomoseq_test_x, tomoseq_test_y, tomoseq_test_z, overwrite = TRUE)

set.seed(116)
n <- 13
p <- 101
k <- 2
q <- 3
is_null       <- rep(FALSE, length = p)
is_null[1:57] <- TRUE

X <- matrix(stats::rnorm(n * q), nrow = n)
B <- matrix(stats::rnorm(q * p), nrow = q)
B[2, is_null] <- 0
Z <- X %*% matrix(stats::rnorm(q * k), nrow = q) +
     matrix(rnorm(n * k), nrow = n)
A <- matrix(stats::rnorm(k * p), nrow = k)
E <- matrix(stats::rnorm(n * p, sd = 1 / 2), nrow = n)
Y <- X %*% B + Z %*% A + E
mout <- mouthwash(Y = Y, X = X, k = k, cov_of_interest = 2,
                  include_intercept = FALSE)

data = rnaseq_data$body
mout <- mouthwash(Y = t(data$l2c), X = data$mod1, k = 1, cov_of_interest = 2,
                  include_intercept = FALSE)
names(mout)
mout$Zhat
residuals = getBatchResiduals(x = data$l2c, tissue = data$tissue, label = "l2c-mouthwash-1d",
                              design = data$design, mod0 = data$mod0, sva = mout$Zhat, col = pull(data$covariates, treatment))

data = rnaseq_data$body
mout2 <- mouthwash(Y = t(data$l2c), X = data$mod1, k = 2, cov_of_interest = 2,
                  include_intercept = FALSE)
mout2$Zhat
residuals = getBatchResiduals(x = data$l2c, tissue = data$tissue, label = "l2c-mouthwash-2d",
                              design = data$design, mod0 = data$mod0, sva = mout2$Zhat, col = pull(data$covariates, treatment))

table <- 
     list_rbind(list(head = mutate(rnaseq_data$head$mwash$l2c$result,
                                   Geneid = rownames(rnaseq_data$head$mwash$l2c$result)),
                    body = mutate(rnaseq_data$body$mwash$l2c$result,
                                   Geneid = rownames(rnaseq_data$body$mwash$l2c$result))), 
               names_to = "tissue") %>% 
     as_tibble %>% 
     relocate(Geneid)

table = left_join(table,
        bitr(table$Geneid, 
             fromType="FLYBASE", 
             toType = "SYMBOL",
             OrgDb = org.Dm.eg.db, 
             drop = FALSE),
        by = c(Geneid = "FLYBASE"))  %>% relocate(Geneid, SYMBOL)
write_csv(table, affix_date("output/HS-ctrl-DE_vicar-head2-body2-table.csv"))




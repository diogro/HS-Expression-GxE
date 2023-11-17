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
write_csv(table, affix_date("output/HS-ctrl-DE_vicar-table.csv"))

data = rnaseq_data$body
mout2 <- mouthwash(Y = t(data$l2c), X = data$mod1, k = 2, cov_of_interest = 2,
                   include_intercept = FALSE)

table2 <- mutate(mout2$result,
                Geneid = rownames(rnaseq_data$body$mwash$l2c$result),
                tissue = "body") %>%
     as_tibble 
table2 = left_join(table2,
        bitr(table$Geneid, 
             fromType="FLYBASE", 
             toType = "SYMBOL",
             OrgDb = org.Dm.eg.db, 
             drop = FALSE),
        by = c(Geneid = "FLYBASE"))  %>% relocate(Geneid, SYMBOL, tissue)
write_csv(table2, affix_date("output/HS-ctrl-DE_vicar_body-2ruv-table.csv"))


mout2$result == 
rnaseq_data$body$mwash$l2c$result

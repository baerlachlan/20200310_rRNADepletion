library(tidyverse)
tibble(
    position = 1:5,
    inKmer = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    s1 = rnorm(5),
    s2 = rnorm(5),
    s3 = rnorm(5)
) %>%
    pivot_longer(cols = starts_with("s"), names_to = "sample", values_to = "coverage") %>%
    mutate(
        riboPerc  = case_when(
            sample == "s1" ~ 5,
            sample == "s2" ~ 8,
            sample == "s3" ~ 10
        )
    ) %>%
    with(
        lm(coverage ~ inKmer*riboPerc)
    ) %>%
    summary()

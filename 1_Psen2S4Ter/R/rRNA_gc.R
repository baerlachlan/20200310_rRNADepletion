library(Biostrings)
library(tidyverse)

rrna <- readDNAStringSet(here("ext_data/rRNA.fa"))
fs <- 300

Views(rrna[[1]], start = seq(1, 42999 - fs, by = fs / 10), width = fs) %>% 
    letterFrequency(letters = "GC", as.prob = TRUE) %>% 
    as_tibble() %>%
    mutate(window = ((seq_along(`G|C`) - 1) * fs / 10) + 1) %>%
    ggplot(aes(window, `G|C`)) + 
    geom_point() +
    geom_smooth(se = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Window start position (bp)", y = "GC content") +
    theme(text = element_text(size = 15))

ggsave("~/phd/thesis/tex/fig/rRNA/rrna_gc.pdf", height = 4, width = 8)

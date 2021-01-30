### Bin approach

#### Dataset 1

```{r binTable_1}
nBins <- list(length = 20, gc = 20)
binTable_1 <- geneGcLen %>%
    left_join(
        pca_1$rotation %>%
            as.data.frame() %>%
            rownames_to_column("gene_id")
    ) %>%
    na.omit() %>%
    as_tibble() %>%
    mutate(
        lengthBins = cut(
            log(aveLen), 
            breaks = quantile(
                log(aveLen), seq(0, nBins$length)/nBins$length
            ),
            labels = paste0("L", seq_len(nBins$length)), 
            include.lowest = TRUE
        ),
        gcBins = cut(
            aveGc, 
            breaks = quantile(
                aveGc, seq(0, nBins$gc) / nBins$gc
            ),
            labels = paste0("GC", seq_len(nBins$gc)), 
            include.lowest = TRUE
        ),
        bothBins = paste(lengthBins, gcBins, sep = "_"),
        bothBins = as.factor(bothBins)
    ) %>%
    group_by(bothBins) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    dplyr::filter(n > 1) %>%
    as_tibble()
```

```{r binFuncs_1}
minLen <- function(x){
    binTable_1 %>% 
        dplyr::filter(lengthBins == x) %>%
        dplyr::select(aveLen) %>% 
        min()
}
maxLen <- function(x){
    binTable_1 %>% 
        dplyr::filter(lengthBins == x) %>%
        dplyr::select(aveLen) %>% 
        max()
}
minGC <- function(x){
    binTable_1 %>% 
        dplyr::filter(gcBins == x) %>%
        dplyr::select(aveGc) %>% 
        min()
}
maxGC <- function(x){
    binTable_1 %>% 
        dplyr::filter(gcBins == x) %>%
        dplyr::select(aveGc) %>% 
        max()
}
```

```{r binRanges_1}
binRanges_1 <- tibble(
    Bin = c(1:20),
    minLen = c(minLen("L1"), minLen("L2"), minLen("L3"), minLen("L4"), 
               minLen("L5"), minLen("L6"), minLen("L7"), minLen("L8"),
               minLen("L9"), minLen("L10"), minLen("L11"), minLen("L12"), 
               minLen("L13"), minLen("L14"), minLen("L15"), minLen("L16"), 
               minLen("L17"), minLen("L18"), minLen("L19"), minLen("L20")),
    maxLen = c(maxLen("L1"), maxLen("L2"), maxLen("L3"), maxLen("L4"), 
               maxLen("L5"), maxLen("L6"), maxLen("L7"), maxLen("L8"),
               maxLen("L9"), maxLen("L10"), maxLen("L11"), maxLen("L12"), 
               maxLen("L13"), maxLen("L14"), maxLen("L15"), maxLen("L16"), 
               maxLen("L17"), maxLen("L18"), maxLen("L19"), maxLen("L20")),
    minGC = c(minGC("GC1"), minGC("GC2"), minGC("GC3"), minGC("GC4"), 
              minGC("GC5"), minGC("GC6"), minGC("GC7"), minGC("GC8"),
              minGC("GC9"), minGC("GC10"), minGC("GC11"), minGC("GC12"), 
              minGC("GC13"), minGC("GC14"), minGC("GC15"), minGC("GC16"), 
              minGC("GC17"), minGC("GC18"), minGC("GC19"), minGC("GC20")),
    maxGC = c(maxGC("GC1"), maxGC("GC2"), maxGC("GC3"), maxGC("GC4"), 
              maxGC("GC5"), maxGC("GC6"), maxGC("GC7"), maxGC("GC8"),
              maxGC("GC9"), maxGC("GC10"), maxGC("GC11"), maxGC("GC12"), 
              maxGC("GC13"), maxGC("GC14"), maxGC("GC15"), maxGC("GC16"), 
              maxGC("GC17"), maxGC("GC18"), maxGC("GC19"), maxGC("GC20"))
) %>%
    mutate(
        minGC = formatC(100*minGC, format = "f", digits = 1),
        maxGC = formatC(100*maxGC, format = "f", digits = 1),
        minLen = formatC(minLen / 1000, format = "f", digits = 1),
        maxLen = formatC(maxLen / 1000, format = "f", digits = 1)
    ) %>%
    unite("Length", minLen, maxLen, sep = "-") %>%
    unite("GC_Content", minGC, maxGC, sep = "-")
```

```{r lmPC1_1}
lmPC1_1 <- binTable_1 %>%  
    lm(PC1 ~ 0 + bothBins, data = .)
```

```{r fdrAdjust_1}
binFDR_PC1 <- summary(lmPC1_1)$coef %>%
    as.data.frame() %>% 
    rownames_to_column("Bin") %>% 
    as_tibble() %>%
    mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) %>%
    mutate(Bin = str_remove(Bin, "bothBins")) %>% 
    separate(Bin, into = c("length", "GC")) %>%
    dplyr::mutate(
        sig = FDR < 0.05,
        p = `Pr(>|t|)`
    ) %>%
    dplyr::select(length, GC, sig, FDR, p)
```

```{r plotBins_1}
coef(lmPC1_1) %>% 
    enframe() %>%
    mutate(name = str_remove(name, "bothBins")) %>% 
    separate(name, into = c("length", "GC")) %>%
    left_join(binFDR_PC1) %>%
    dplyr::filter(sig) %>%
    mutate(
        length = str_remove(length, "L") %>% as.integer(),
        GC = str_remove(GC, "GC") %>% as.integer()
    ) %>%
    ggplot(aes(length, GC)) +
    geom_point(aes(colour = value, size = -log10(p))) +
    scale_colour_gradient2(
        low = "red", mid = "white", high = "blue", midpoint = 0
    ) +
    scale_x_continuous(
        breaks = seq_len(nBins$length),
        labels = binRanges_1$Length
    ) +
    scale_y_continuous(
        breaks = seq_len(nBins$gc),
        labels = binRanges_1$GC_Content
    ) +
    scale_size_continuous(guide = FALSE) +
    labs(
        x = "Gene length (kb)",
        y = "GC content (%)",
        colour = paste("Contribution", "to PC1", sep = "\n"),
        size = "Significance",
        title = "E-GEOD-71609"
    ) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = -45, hjust = 0.1, vjust = 0.5
        )
    )
```

#### Dataset 2

```{r binTable_2}
nBins <- list(length = 20, gc = 20)
binTable_2 <- geneGcLen %>%
    left_join(
        pca_2$rotation %>%
            as.data.frame() %>%
            rownames_to_column("gene_id")
    ) %>%
    na.omit() %>%
    as_tibble() %>%
    mutate(
        lengthBins = cut(
            log(aveLen), 
            breaks = quantile(
                log(aveLen), seq(0, nBins$length)/nBins$length
            ),
            labels = paste0("L", seq_len(nBins$length)), 
            include.lowest = TRUE
        ),
        gcBins = cut(
            aveGc, 
            breaks = quantile(
                aveGc, seq(0, nBins$gc) / nBins$gc
            ),
            labels = paste0("GC", seq_len(nBins$gc)), 
            include.lowest = TRUE
        ),
        bothBins = paste(lengthBins, gcBins, sep = "_"),
        bothBins = as.factor(bothBins)
    ) %>%
    group_by(bothBins) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    dplyr::filter(n > 1) %>%
    as_tibble()
```

```{r binFuncs_2}
minLen <- function(x){
    binTable_2 %>% 
        dplyr::filter(lengthBins == x) %>%
        dplyr::select(aveLen) %>% 
        min()
}
maxLen <- function(x){
    binTable_2 %>% 
        dplyr::filter(lengthBins == x) %>%
        dplyr::select(aveLen) %>% 
        max()
}
minGC <- function(x){
    binTable_2 %>% 
        dplyr::filter(gcBins == x) %>%
        dplyr::select(aveGc) %>% 
        min()
}
maxGC <- function(x){
    binTable_2 %>% 
        dplyr::filter(gcBins == x) %>%
        dplyr::select(aveGc) %>% 
        max()
}
```

```{r binRanges_2}
binRanges_2 <- tibble(
    Bin = c(1:20),
    minLen = c(minLen("L1"), minLen("L2"), minLen("L3"), minLen("L4"), 
               minLen("L5"), minLen("L6"), minLen("L7"), minLen("L8"),
               minLen("L9"), minLen("L10"), minLen("L11"), minLen("L12"), 
               minLen("L13"), minLen("L14"), minLen("L15"), minLen("L16"), 
               minLen("L17"), minLen("L18"), minLen("L19"), minLen("L20")),
    maxLen = c(maxLen("L1"), maxLen("L2"), maxLen("L3"), maxLen("L4"), 
               maxLen("L5"), maxLen("L6"), maxLen("L7"), maxLen("L8"),
               maxLen("L9"), maxLen("L10"), maxLen("L11"), maxLen("L12"), 
               maxLen("L13"), maxLen("L14"), maxLen("L15"), maxLen("L16"), 
               maxLen("L17"), maxLen("L18"), maxLen("L19"), maxLen("L20")),
    minGC = c(minGC("GC1"), minGC("GC2"), minGC("GC3"), minGC("GC4"), 
              minGC("GC5"), minGC("GC6"), minGC("GC7"), minGC("GC8"),
              minGC("GC9"), minGC("GC10"), minGC("GC11"), minGC("GC12"), 
              minGC("GC13"), minGC("GC14"), minGC("GC15"), minGC("GC16"), 
              minGC("GC17"), minGC("GC18"), minGC("GC19"), minGC("GC20")),
    maxGC = c(maxGC("GC1"), maxGC("GC2"), maxGC("GC3"), maxGC("GC4"), 
              maxGC("GC5"), maxGC("GC6"), maxGC("GC7"), maxGC("GC8"),
              maxGC("GC9"), maxGC("GC10"), maxGC("GC11"), maxGC("GC12"), 
              maxGC("GC13"), maxGC("GC14"), maxGC("GC15"), maxGC("GC16"), 
              maxGC("GC17"), maxGC("GC18"), maxGC("GC19"), maxGC("GC20"))
) %>%
    mutate(
        minGC = formatC(100*minGC, format = "f", digits = 1),
        maxGC = formatC(100*maxGC, format = "f", digits = 1),
        minLen = formatC(minLen / 1000, format = "f", digits = 1),
        maxLen = formatC(maxLen / 1000, format = "f", digits = 1)
    ) %>%
    unite("Length", minLen, maxLen, sep = "-") %>%
    unite("GC_Content", minGC, maxGC, sep = "-")
```

```{r lmPC1_2}
lmPC1_2 <- binTable_2 %>%  
    lm(PC1 ~ 0 + bothBins, data = .)
```

```{r fdrAdjust_2}
binFDR_PC1 <- summary(lmPC1_2)$coef %>%
    as.data.frame() %>% 
    rownames_to_column("Bin") %>% 
    as_tibble() %>%
    mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) %>%
    mutate(Bin = str_remove(Bin, "bothBins")) %>% 
    separate(Bin, into = c("length", "GC")) %>%
    dplyr::mutate(
        sig = FDR < 0.05,
        p = `Pr(>|t|)`
    ) %>%
    dplyr::select(length, GC, sig, FDR, p)
```

```{r plotBins_2}
coef(lmPC1_2) %>% 
    enframe() %>%
    mutate(name = str_remove(name, "bothBins")) %>% 
    separate(name, into = c("length", "GC")) %>%
    left_join(binFDR_PC1) %>%
    dplyr::filter(sig) %>%
    mutate(
        length = str_remove(length, "L") %>% as.integer(),
        GC = str_remove(GC, "GC") %>% as.integer()
    ) %>%
    ggplot(aes(length, GC)) +
    geom_point(aes(colour = value, size = -log10(p))) +
    scale_colour_gradient2(
        low = "red", mid = "white", high = "blue", midpoint = 0
    ) +
    scale_x_continuous(
        breaks = seq_len(nBins$length),
        labels = binRanges_2$Length
    ) +
    scale_y_continuous(
        breaks = seq_len(nBins$gc),
        labels = binRanges_2$GC_Content
    ) +
    scale_size_continuous(guide = FALSE) +
    labs(
        x = "Gene length (kb)",
        y = "GC content (%)",
        colour = paste("Contribution", "to PC1", sep = "\n"),
        size = "Significance",
        title = "E-GEOD-72322"
    ) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = -45, hjust = 0.1, vjust = 0.5
        )
    )
```
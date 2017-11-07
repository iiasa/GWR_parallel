# GWR_parallel


Time using 32 cores. Tine in secs.

```{r}
tibble::tribble(    ~N,      ~Time,     ~Cores,
                   100,      8.480, "par_ggwr",
                   500,     18.649, "par_ggwr",
                  1000,     38.281, "par_ggwr",
                  2000,     84.007, "par_ggwr",
                  5000,    213.155, "par_ggwr",
                   100,      6.831,     "ggwr",
                   500,     27.211,     "ggwr",
                  1000,     53.249,     "ggwr",
                  2000,    101.376,     "ggwr",
                  5000,    259.116,     "ggwr") %>% 
  ggplot2::ggplot(aes(x = N, y = Time, group = Cores, colour = Cores)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm")
```
# GWR_parallel


Time using 32 cores. Tine in secs.

```{r}
tibble::tribble(    ~N,      ~Time,                  ~Version,
                   100,      8.480,                "par_ggwr",
                   500,     18.649,                "par_ggwr",
                  1000,     38.281,                "par_ggwr",
                  2000,     84.007,                "par_ggwr",
                  5000,    213.155,                "par_ggwr",
                   100,      6.831,                    "ggwr",
                   500,     27.211,                    "ggwr",
                  1000,     53.249,                    "ggwr",
                  2000,    101.376,                    "ggwr",
                  5000,    259.116,                    "ggwr",
                   100,      3.871, "par_ggwr_max_dist_4000km",
                   500,      5.668, "par_ggwr_max_dist_4000km",
                  1000,      7.700, "par_ggwr_max_dist_4000km",
                  2000,     17.550, "par_ggwr_max_dist_4000km",
                  5000,     40.782, "par_ggwr_max_dist_4000km") %>% 
  ggplot2::ggplot(aes(x = N, y = Time, group = Version, colour = Version)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm")
```
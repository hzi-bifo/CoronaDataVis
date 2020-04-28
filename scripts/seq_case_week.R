library(tidyverse)
library(zoo)
library(ggsci)
library(grid)

country_rename_map <- c(`Democratic Republic of the Congo`='Congo',
                        `"Korea, South"`= 'South Korea',
                        `Korea, South`= 'South Korea',
                        `Congo (Brazzaville)`= 'Congo',
                        `Congo (Kinshasa)`= 'Congo',
                        `Czech Republic`='Czechia',
                        `Hong Kong`='China',
                        `Taiwan*`='China',
                        Taiwan='China',
                        US='USA')

seq_meta <- read_tsv("../data/metadata_fullgenome_filtered.tsv") %>%
  mutate(country=ifelse(country %in% names(country_rename_map), 
                                 country_rename_map[country],
                                 country))
case_meta <- read_csv("../data/time_series_covid19_confirmed_global.csv") %>% 
  mutate(`Country/Region`=ifelse(`Country/Region` %in% names(country_rename_map), 
                                 country_rename_map[`Country/Region`],
                                 `Country/Region`))

seq_meta %>% 
  group_by(country, date) %>%
  tally() %>%
  arrange(date) ->
  seq_country_date

case_meta[,-c(1,3,4)] %>% 
  group_by(`Country/Region`) %>%
  summarise_all(sum) %>%
  gather(date, case, -`Country/Region`) %>%
  mutate(date=as.Date(date, "%m/%d/%y")) %>%
  arrange(date) ->
  case_country_date

colnames(case_country_date) <- c("country", "date", "case")


seq_country_date %>% 
  group_by(country) %>% 
  summarise(total_seqs=sum(n)) -> final_seq_country

case_country_date %>% group_by(country) %>% summarise(total_cases=last(case)) ->
  final_case_country

left_join(final_seq_country, final_case_country, by=c("country")) %>% 
  mutate(fraction=total_seqs/total_cases) %>%
  arrange(desc(fraction)) -> final_seq_case

final_seq_case$country <- factor(final_seq_case$country, levels = final_seq_case$country)


barplot_rank <- ggplot(final_seq_case, aes(country, fraction)) + geom_bar(stat = 'identity') + 
  theme_bw(base_size = 20) +
  xlab("") +
  ylab("# Sequences / # Cases") +
  theme(axis.text.x=element_text(angle = 70, hjust=1), text = element_text(size = 25))

ggsave("../figures/seq_case_fraction_rank.pdf", plot = barplot_rank, width = 17)


week_seq_case_frac <- function(seq, case){

  total <- tail(seq, 1)$n
  aligned_head_seq <- bind_rows(tibble(country=seq$country[1], date=case$date[1], n=0), seq)
    
   
  aligned_tail_seq <- bind_rows(aligned_head_seq, tibble(country=seq$country[1], date=tail(case, 1)$date, n=total))
  
  aligned_tail_seq$week <- format(aligned_tail_seq$date, "%Y-W%U")
  case$week <- format(case$date, "%Y-W%U")
  
  agg_week_seq <- aligned_tail_seq %>% group_by(country, week) %>% summarise(week_n=last(n))
  agg_week_case <- case %>% group_by(country, week) %>% summarise(week_case=last(case))
  
  agg_week_case_seq <- full_join(agg_week_case, agg_week_seq, by=c("country", "week")) %>% 
    do(na.locf(.)) %>% mutate(frac=ifelse(week_case==0, 0, week_n/week_case)) 
  
  agg_week_case_seq[agg_week_case_seq$frac > 1, 'frac'] <- 1
  
  return(agg_week_case_seq)
}

country_list <- c("Germany", "USA", 'China', 
                  'Iceland', 'France',
                  'Italy', 'Netherlands',  
                  'Luxembourg',  
                  'United Kingdom')#, 'Spain')

create_case_seq_df <- function(seq, case, country_list)
{
  seq_case_list <- list()
  for (country_name in country_list){
    seq %>% dplyr::filter(country==country_name) %>% 
      mutate(n=cumsum(n)) -> country_seq
    case %>% dplyr::filter(country==country_name) -> country_case
    
    country_seq_case <- week_seq_case_frac(country_seq, country_case) 
    seq_case_list[[country_name]] <- country_seq_case
  }
  
  return(do.call(rbind, seq_case_list))

}

case_seq_df <- create_case_seq_df(seq_country_date, case_country_date, country_list) %>% 
  filter(week >= "2020-W08")



p <- ggplot(case_seq_df, aes(week, frac, group=country)) + 
  geom_line(aes(linetype=country, color=country), size=2) +
  geom_point(aes(color=country), size=5) +
  geom_text(data = subset(final_seq_case, country %in% country_list), 
            aes(label = total_seqs, colour = country, x = Inf, y = fraction), hjust = -.1,
            size = 7) +

  theme_bw(base_size = 27) +
  theme(axis.text.x=element_text(angle = 70, hjust=1), 
        text = element_text(size = 27), 
        legend.position="right",
        legend.justification="right",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(10,10,10,35)) +
  ylab("# Sequences / # Cases") + 
  scale_color_ucscgb() +
  scale_y_log10()

gt <- ggplotGrob(p)
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

ggsave("../figures/seq_case_fraction_curve.pdf", plot = gt, width = 14, height = 8)




#' ---
#' title: "Class 05 Data Visualization"
#' author: "Monica Lin (PID: A15524235)"
#' ---

# Anything you type in this file must be run in the console below
# The code typed here will not auto-run

# Lets start with a scatterplot
cars

# We have already installed ggplot, but if not, type the following line:
# install.packages("ggplot2")
# Before we can use it, we need to load it up!
library(ggplot2)

# Every ggplot has a data + aes + geoms layer
ggplot(data=cars) + 
  aes(x=speed, y=dist) +
  geom_point() + 
  geom_smooth()

# Change to a linear model
p <- ggplot(data=cars) + 
  aes(x=speed, y=dist) + 
  geom_point() + 
  geom_smooth(method="lm")

p + labs(title="My nice plot", 
         x="Speed (MPH)", 
         y="Stopping Distance (ft)",
         caption="Dataset: 'cars'")

# Base graphics is shorter
plot(cars)

# Can add B/W theme
ggplot(cars) + 
  aes(x=speed, y=dist) +
  geom_point() +
  labs(title="Speed and Stopping Distances of Cars",
       x="Speed (MPH)", 
       y="Stopping Distance (ft)",
       subtitle = "Monica Lin",
       caption="Dataset: 'cars'") +
  geom_smooth(method="lm", se=FALSE) +
  theme_bw()

# Adding more plot aesthetics through aes()
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

# How many genes in this dataset?
nrow(genes)

# Use colnames() and ncol() on the genes data.frame
colnames(genes)
ncol(genes)

# Use table() on the State column of this data.frame
# How many 'up' regulated genes are there?
table(genes$State)

# What % are up/down regulated?
table(genes$State) / nrow(genes) * 100
# Use the round() function to auto display sig figs
round(table(genes$State) / nrow(genes) * 100, 2)
# Alternative way to round
prec <- table(genes$State) / nrow(genes) * 100
round(prec, 3)

# Display the genes data
ggplot(genes) +
  aes(x=Condition1, y=Condition2) +
  geom_point()

# State column tells stat sig
p <- ggplot(genes) + 
  aes(x=Condition1, y=Condition2, col=State) + 
  geom_point()
p

# Change to user-friendly colors
p + scale_color_manual(values=c("blue", "gray", "red"))

# Add labels
p + scale_color_manual(values=c("blue", "gray", "red")) +
  labs(title="Gene Expression Changes Upon Drug Treatment",
       x="Control (no drug)",
       y="Drug Treatment")

# OPTIONAL: Going Further
# gapminder dataset contains econ + demographic data on various countries since 1952
url <-  "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"
gapminder <- read.delim(url)

# Plot animation!
# Install gganimate and gifski to make animated gifs
install.packages("gifski")
install.packages("gganimate")

# Install and call up gapminder data, and gganimate
install.packages("gapminder")
library(gapminder)
library(gganimate)

#Setup nice regular ggplot of the gapminder data
ggplot(gapminder) + 
  aes(gdpPercap, lifeExp, size=pop, color=country) +
  geom_point(alpha=0.7, show.legend=FALSE) +
  scale_color_manual(values=country_colors) +
  scale_size(range=c(2,12)) +
  scale_x_log10() +
  facet_wrap(~continent) + 
  labs(title="Year:{frame_time}", x="GDP per capita", y="life expectancy") +
  transition_time(year) +
  shadow_wake(wake_length=0.1, alpha=FALSE)

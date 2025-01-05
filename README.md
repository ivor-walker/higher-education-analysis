# Introduction
A lecturer at St Andrews teaches a ten-week intensive statistics course to university students. She is interested to know whether students who took Higher mathematics improve their performance at a different rate than those without. She has collected data from 80 students, 40 with Higher mathematics and 40 without, selected at random from the cohort from each country in last year’s class, and recorded their improvement on each assessment each week. 

# Model
After discussion, we decided to fit the following model:

yit = μ + β exp(−(δ + γsi)t) + ϵit

where:
- yit is the improvement of student i on assessment t,
- μ is the overall improvement, 
- β is the overall rate of improvement, 
- δ is the rate of decay of improvement, 
- γ is the difference in rate of improvement between students with and without Higher mathematics
- si is an indicator variable for whether student i took Higher mathematics.- ϵit is an error term assumed to be normally distributed with mean 0 and variance σ2.

# Solution
I obtained starting values for each of these parameters from the data and model (see R comments). Then, I used the Gauss-Newton algorithm to produce estimates for the parameters such that they minimise the sum of squared errors between the observed and predicted values of the dependent variable. By finding the optimal values of γ, I concluded that students with Higher mathematics improved at a different rate than those without. I was able to provide other insights into her course, such as how the rate of improvement changes over time and the overall improvement of her students.

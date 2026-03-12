# Stochastic processes

## Arithmetic Brownian Motion 
$$dX_t = \mu dt + \sigma dW_t$$
For a time step $\Delta t_i = (s-t), the exact dicretization is given by: 
$$X_s = X_t + \mu\Delta_t+ \sigma \sqrt{\Delta_t}Z_t \Longleftrightarrow \Delta X_i = \mu\Delta t_i + \sigma \sqrt{\Delta t_i }Z_t$$
With $Z_t \sim N(0,1)$. Therefore transition probability is given by: 
$$\Delta X_i  \sim N(\mu\Delta t_i, \sigma^2\Delta t_i)$$
The likelihood is therefore given by: 
$$ L = \prod_{i=1}^n \frac{1}{\sqrt{2\pi \sigma^2\Delta t_i}} \exp{\left[-\frac{(\Delta X_i - \mu \Delta t_i)^2}{2\sigma^2 \Delta t_i}\right]}$$
It follows that: 
$$\log L = -\frac{1}{2} \sum_{i=1}^n \log(2\pi\sigma^2\Delta t_i) + \frac{(\Delta X_i - \mu \Delta t_i)^2}{\sigma^2 \Delta t_i}$$

$$\partial_\mu \log L = 0 \Longleftrightarrow \hat{\mu} = \frac{\sum_{i=1}^n \Delta X_i}{\sum_{i=1}^n \Delta t_i}$$
$$\partial_\sigma \log L = 0 \Longleftrightarrow \hat{\sigma}^2 = \frac{1}{n} \sum_{i=1}^n \frac{(\Delta X_i -\mu\Delta t_i)^2}{\Delta t_i}$$
Additionnaly, by transformation of the discretization, we can see that: 
$$\frac{\Delta X_i}{\sqrt{\Delta t_i}} = \mu \sqrt{\Delta t_i} + \sigma Z_t\$$
Which becomes an OLS-type problem. 
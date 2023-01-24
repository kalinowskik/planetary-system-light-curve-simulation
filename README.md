# planetary system light curve simulation
 code that allows to generate synthetic light curves of planetary systems, could be used to study detectability of planets through fourier transform of the light curve


## Model

### The observed system

The system is modelled following Charpinet et al 2011, Supplementary Information, section D (page 18) (https://www.nature.com/articles/nature10631#Sec2).

A star with temperature
$T_*$
emits black body radiation with power $ L_* $. The orbit of each planet in the system is circular and has a radius of $d_i$. The radius of the $i$-th planet is $r_i$. A portion of the star's radiation equal to the albedo $α_i$ of the planet falling on it is immediately reflected by its surface. The remaining part of the falling radiation is absorbed and heats the planet. All absorbed energy is re-emitted. The bright side of the planet has a temperature $T_{ib}$. It is determined from the formula
$$T_{ib} = (1 - α_i) / (8πσ) * L_* / d_i^2 * 1 / (1 + β_i^4)^{1/4},$$
where $β$ is a constant defined to describe the degree of energy redistribution from the star to the planet's surface.

The dark side has a temperature $T_{id} = βT_{ib}$.

The power of radiation received by the detector depends on the transmission profile of the filter $B(λ)$ used. If the radiation comes from a perfect black body, its flux can be expressed by the formula.

$$
    F(T)=\int_0^{\infty}B_{\lambda}(T)L_{\lambda}d\lambda\text{,}
$$

where $λ$ is the wavelength, $B_λ$ is the Planck distribution for a body at temperature $T$, and $L$ is the transmission curve of the filter used.

The power of the planet on the side facing the observer can be described by the formula:

$$L_i = L_* α_i / 8 * r_i^2 / d_i^2 + 2π r_i^2 (F(T_b) + F(T_d)) + (L_* α_i / 8 * r_i^2 / d_i^2 + 2π r_i^2 (F(T_b) - F(T_d))) \sin(i_i) \sin(ϕ(t))$$
where $i$ is the inclination of the planet's orbit and $ϕ(t)$ is the anomaly dependent on time.

### Observation Mechanics
The observations consisted of an immediate measurement of the total brightness of the planets (i.e. substitution of data into formula for $L_i$) at certain intervals of time, i.e. the sampling period, $τ_p$. The observations were performed for the duration of the campaign, $τ_{obs}$.

During actual observational campaigns, astronomers are faced with random events that make measurements impossible, such as cloud cover. As a result, it is impossible to make some of the planned observations. This was reflected in the simulation: random observations were selected without replacement and removed. The parameter $κ$ describes the portion of all observations that were removed. For example, if $κ = 0$, all observations were saved, and if $κ = 1$, no measurements were made.

In addition, light pollution from the sun is a problem for ground-based observers, making observations impossible from dawn to dusk. This was reflected in the simulation by defining the length of the day, $P_d$, and the length of the night, $P_n$. It was checked for each measurement whether it fell during the night or during the day, and in the latter case, the measurement was removed. Therefore, if $P_d = P_n$, no measurements were removed.

Additionally, the simulation took into account typical noise for physical measurements, i.e. a value drawn from a normal distribution with a standard deviation of $σ$ was added to each measurement value.

It was assumed that the phase of all planets was the same at the start of the observation.
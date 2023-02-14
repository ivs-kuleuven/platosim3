import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from platosim.lightcurve import LightCurve
from transitleastsquares import transitleastsquares


lc = LightCurve('test.ftr')
df = lc.data()

# Cut data
x = 50000000
time, flux = df.time.iloc[:x], df.flux.iloc[:x]

# Run TLS
model = transitleastsquares(df.time, df.flux)
results = model.power(R_star_min=0.9, R_star_max=1.1,
                      M_star_min=0.9, M_star_max=1.1,
                      period_min=350., period_max=380.,
                      use_threads=10)


print('Period', format(results.period, '.5f'), 'd')
print(len(results.transit_times), 'transit times in time series:', \
        ['{0:0.5f}'.format(i) for i in results.transit_times])
print('Transit depth', format(results.depth, '.5f'))
print('Best duration (days)', format(results.duration, '.5f'))
print('Signal detection efficiency (SDE):', results.SDE)



plt.figure()
ax = plt.gca()
ax.axvline(results.period, alpha=0.4, lw=3)
plt.xlim(numpy.min(results.periods), numpy.max(results.periods))
for n in range(2, 10):
    ax.axvline(n*results.period, alpha=0.4, lw=1, linestyle="dashed")
    ax.axvline(results.period / n, alpha=0.4, lw=1, linestyle="dashed")
plt.ylabel(r'SDE')
plt.xlabel('Period (days)')
plt.plot(results.periods, results.power, color='black', lw=0.5)
plt.xlim(0, max(results.periods))

plt.figure()
plt.plot(results.model_folded_phase, results.model_folded_model, color='red')
plt.scatter(results.folded_phase, results.folded_y, color='blue', s=10, alpha=0.5, zorder=2)
plt.xlim(0.48, 0.52)
plt.ticklabel_format(useOffset=False)
plt.xlabel('Phase')
plt.ylabel('Relative flux');


bins = 500
bin_means, bin_edges, binnumber = stats.binned_statistic(
    results.folded_phase,
    results.folded_y,
    statistic='mean',
    bins=bins)
bin_stds, _, _ = stats.binned_statistic(
    results.folded_phase,
    results.folded_y,
    statistic='std',
    bins=bins)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2

plt.plot(results.model_folded_phase, results.model_folded_model, color='red')
plt.scatter(results.folded_phase, results.folded_y, color='blue', s=10, alpha=0.5, zorder=2)
plt.errorbar(
    bin_centers,
    bin_means,
    yerr=bin_stds/2,
    xerr=bin_width/2,
    marker='o',
    markersize=8,
    color='black',
    #capsize=10,
    linestyle='none')
plt.xlim(0.48, 0.52)
plt.ticklabel_format(useOffset=False)
plt.xlabel('Phase')
plt.ylabel('Relative flux');

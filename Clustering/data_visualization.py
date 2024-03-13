# %% Data visualize
import matplotlib.pyplot as plt

exec(open("data_processing.py").read())

# Compare healthy and unhealthy
fig, axes = plt.subplots(len(Feature), 2, sharex=True, sharey='row')

axes[0, 0].set_title('Healthy')
axes[0, 1].set_title('Post-MI')

for i, feature in enumerate(Feature):
    time_Ht, feat_Ht = join_stimulus(dataframes_Ht, 1, i)
    axes[i, 0].plot(time_Ht, feat_Ht)
    axes[i, 0].set_ylabel(feature)

    time_MI, feat_MI = join_stimulus(dataframes_MI, 3, i)
    axes[i, 1].plot(time_MI, feat_MI)

fig.supxlabel('Time (s)')
plt.tight_layout()
plt.savefig('Figure_1.png')
plt.show()

# Individual stimulus

fig, axes = plt.subplots(len(Feature), len(Stimulus), sharex='col', sharey='row')

for i, feature in enumerate(Feature):
    axes[i, 0].set_ylabel(feature)

    for j, stim in enumerate(Stimulus):
        time_Ht, feat_Ht = access_df(dataframes_Ht, 1, j, i)
        time_MI, feat_MI = access_df(dataframes_MI, 3, j, i)

        axes[i, j].plot(time_Ht, feat_Ht)
        axes[i, j].plot(time_MI, feat_MI)

        axes[len(Feature) - 1, j].set_xlabel(stim)

fig.supxlabel('Time (s)')

fig.legend(['Healthy', 'Post-MI'], loc='upper center', ncols=2, framealpha=1)

plt.tight_layout()
plt.savefig('Figure_2.png')
plt.show()

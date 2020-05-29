import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

hidden = [
           [[68, 14], [14, 33]],
           [[71, 11], [15, 32]],
           [[77, 5], [20, 27]],
           [[73, 9], [17, 30]]
          ]

fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
cbar_ax = fig.add_axes([1.01, .3, .03, .4])

for i, ax in enumerate(axs.flat):
    sns.heatmap(hidden[i], ax=ax, cbar=i == 0,
                cbar_ax = None if i else cbar_ax, annot=True,
               cmap="Blues")
    ax.tick_params(axis='both', which='both', length=0)
    if i == 0:
        ax.set_xlabel(f'{i+1} hidden layer')
    else:
        ax.set_xlabel(f'{i+1} hidden layers')

fig.tight_layout()
fig.text(0, 0.45, 'True label', rotation='vertical')
fig.text(0.45, 0, 'Predicted label')

plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
plt.show()

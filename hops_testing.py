import matplotlib.pyplot as plt
import numpy as np

from MegaBeer.model import hop_utilization as hu
# hop utilization plot

t_arr = np.linspace(0., 90., 1000)
G_arr = np.linspace(1., 1.1, 1000)
cmap = plt.get_cmap('cool')

t_grid, G_grid, u_grid = hu.TinsethModel.utilization_grid(t_arr, G_arr)
plt.pcolormesh(t_grid, G_grid, u_grid, cmap=cmap, vmin=0., vmax=np.max(u_grid))
plt.xlim([0., 90.])
plt.ylim([1., 1.1])
plt.xlabel('Time (min)')
plt.ylabel('Gravity')
plt.colorbar()
plt.savefig('hop_utilization.png')
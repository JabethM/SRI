import matplotlib.pyplot as plt

# Sample data
x = range(10)
y1 = [i**2 for i in x]
y2 = [i*10 for i in x]

# Create first plot
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('X')
ax1.set_ylabel('Y1', color=color)
ax1.plot(x, y1, color=color)
ax1.tick_params(axis='y', labelcolor=color)

# Create second plot sharing the same x-axis
ax2 = ax1.twinx()
color = 'tab:blue'
ax2.set_ylabel('Y2', color=color)
ax2.plot(x, y2, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.show()
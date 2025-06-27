import matplotlib.pyplot as plt

# Use a default font that supports English
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['axes.unicode_minus'] = False

# Data
people = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
r = [3, 3, 3, 4, 2, 2, 2, 1]

foods = ['Fruit', 'Meat', 'Drink', 'Cake', 'Food']
c = [4, 2, 6, 4, 4]

# Plot
plt.figure(figsize=(10, 5))

# Subplot 1: Consumption per person
plt.subplot(1, 2, 1)
plt.bar(people, r, color='skyblue')
plt.title('Individual Food Consumption')
plt.xlabel('Person')
plt.ylabel('Amount')

# Subplot 2: Available food portions
plt.subplot(1, 2, 2)
plt.bar(foods, c, color='lightcoral')
plt.title('Available Food Portions')
plt.xlabel('Food Type')
plt.ylabel('Portions')

# Save and show
plt.tight_layout()
plt.savefig("party_distribution.png")
plt.show()

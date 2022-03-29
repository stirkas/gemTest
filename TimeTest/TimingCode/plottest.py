import matplotlib.pyplot as plt
  
x = []
for line in open('testfile.txt', 'r'):
    lines = [i for i in line.split()]
    x.append(lines[0])

plt.title("time per input")
plt.plot(x, marker = 'o', c = 'g')
  
plt.show()
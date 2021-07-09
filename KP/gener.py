import numpy as np

N = 5000
print(N)
b = np.random.randint(1,500,size=(N,N))
b_symm = (b + b.T)/2

for i in range(N):
	for j in range(i,N):
		if (i != j):
			rand = np.random.random()
			if rand < 0.8:
				b_symm[i][j] = 0.0
				b_symm[j][i] = 0.0
sum = 0.0
for i in range(N):
	for j in range(N):
		sum += b_symm[j][i]
	b_symm[i][i] = sum *2 + 15
	sum = 0

for i in range(N):
	for j in range(N):
		print(b_symm[i][j], end=' ')
	print()


b = np.random.randint(0,500, size=N)
for q in b:
	print(q, end=' ')
print()
print('0.001')

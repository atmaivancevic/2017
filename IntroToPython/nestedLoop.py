import turtle

# The outer loop creates a square
# The inner loop specifies the lines to be dashed, 
# with the dashed increasing in size

for _ in range(4):
	for i in range(10):
		turtle.forward(5*i)
		turtle.penup()
		turtle.forward(2)
		turtle.pendown()
	turtle.forward(100)
	turtle.left(90)


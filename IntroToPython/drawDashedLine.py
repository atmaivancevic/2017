import turtle

# draw a simple dashed line
#for i in range(10):
#	turtle.forward(15)
#	turtle.penup()
#	turtle.forward(5)
#	turtle.pendown()

#turtle.exitonclick()

# now make the dashed become larger as the line progresses
# basically utilise i to increase the step size by 2i each time
for i in range(10):
	print(i)
	turtle.forward(15+i*2)
	turtle.penup()
	turtle.forward(5)
	turtle.pendown()

turtle.exitonclick()

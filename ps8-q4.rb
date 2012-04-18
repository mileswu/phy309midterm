n = 1000
r0 = 0.0
a = 2.0
k = 2.0

t = 0.0
c = 0.0
t_step = 2.0*Math::PI / n.to_f
x_prev = nil
y_prev = nil

while(t < 2.0*Math::PI) 
	x = r0 + a*Math.cos(t)
	y = k * a*Math.sin(t)

	if x_prev and y_prev
		c += Math.sqrt((x-x_prev)**2 + (y-y_prev)**2)
	end
	x_prev = x
	y_prev = y

	t += t_step;
end

puts c

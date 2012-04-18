a = 2.0
k = 2.0

def mine(a, k, delta = 0.0)
	n = 1000
	t = 0.0
	c = 0.0
	t_step = 2.0*Math::PI / n.to_f
	x_prev = nil
	y_prev = nil

	while(t < 2.0*Math::PI) 
		x = a*Math.cos(t + delta*Math.sin(t))
		y = k * a*Math.sin(t)

		if x_prev and y_prev
			c += Math.sqrt((x-x_prev)**2 + (y-y_prev)**2)
		end
		x_prev = x
		y_prev = y

		t += t_step;
	end
	return c
end

def simple(a, k)
	return 2.0*Math::PI*a*Math.sqrt((1.0 + k*k)/2.0)
end

def ramanujan(a, k)
	return Math::PI*a*(3.0*(1.0+k) - Math.sqrt((1.0 + 3.0*k) * (3.0 + k)))
end

def qb
	k = 1.0
	a = 2.0

	while(k < 2.6)
		puts "#{sprintf("%.1f", k)} & #{sprintf("%.4f", mine(a,k))} & #{sprintf("%.4f",simple(a,k))} & #{sprintf("%.4f",ramanujan(a,k))} \\\\"
		k += 0.3
	end
end

def qc
	delta = 0
	while(delta < 0.81)
		puts "#{sprintf("%.2f", delta)} & #{sprintf("%.4f", mine(2.0, 1.9, delta))} \\\\"
		delta += 0.16
	end
end

qb
puts ""
qc

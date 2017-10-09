import timeit
start_time = timeit.default_timer()

var('x,y,B,Z')
#var('nu_ai,nu_ei,nu_ii,aspect')

#B = (nu_ai * aspect^(3/2) * nu_ei) / (nu_ii * sqrt(x))
y_func = B / ((y + Z/sqrt(x))^2 + B^2)
#assume(aspect > 0)
assume(x > 0)

inner_yintegral = integrate(y_func, y, -1, 1)
#inner_yintegral.simplify_full()

answer_file = file('/home/kabv/Documents/Masters/Graduation_Project/Model_Source/xi_integral_answer.txt', "w")

answer_file.write("Inner y integral is = "+str(inner_yintegral)+"\n")

# New inner y integral after ALGEBRA
assume(B*(x-Z*sqrt(x)) > 0)
inner_integral = B*arctan(2*Z*sqrt(x) / (B + 1)*x - Z^2)

full_integral = inner_yintegral*x^2*exp(-x)
expanded_integral = maxima(inner_yintegral).powerseries(x, 0)._sage_()
print expanded_integral
#
answer_file.write("Full integral: "+str(expanded_integral)+"\n")

elapsed_time = timeit.default_timer() - start_time
answer_file.write("It took "+str(elapsed_time)+" seconds\n")

answer_file.close()


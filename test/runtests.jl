using RationalExtensions
using Base.Test

# write your own tests here
@test Rad(0,0,1) == 0
@test (1 + 2*Sqrt(3))//(2+3*Sqrt(3)) == 16//23 - Sqrt(3)//23
@test Rad(-big(3)//4,1//2,5) == -3//4 + Sqrt(5)//2
@test Rad(big(1),2,3) == Rad(1,2,3)
@test (r = Sqrt(200); (r.a,r.b,r.n) == (0,10,2))
@test 3//4 * Sqrt(3)//4 + Sqrt(3) == 19Sqrt(3)//16
@test true*(1-2Sqrt(5)) == 1 - 2Sqrt(5)
@test Sqrt(3)/Sqrt(5) ≈ 0.7745966692414833
@test isa(Sqrt(3)/Sqrt(5),Float64)
@test 2/Sqrt(5) ≈ 0.8944271909999159
@test -(1+2Sqrt(3)) == -1 - 2Sqrt(3)
@test big(3//4) + 2//3*Sqrt(1) == 17//12
@test convert(Rad,3//4) == 3//4
@test isa(convert(Rad,3//4),Rad{Rational{Int64},Int64})
@test repr(1-Sqrt(5)//2) == "1 - √(5)//2"

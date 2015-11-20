function xxx(a::Array{Float64,1})
    a=a+1
end

b=[1.0,2.0,3.0]

xxx(b)

b="abc"
symbol(b)

eval(Expr(:(=),symbol(b),123))
abc

type hello
  a::Float64
  function hello(m,n)
    a=1.5
    new(a)
  end
end

go=hello(0.5)
go.a

Expr(:(=),symbol(b),123)

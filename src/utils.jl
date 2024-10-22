# Originally contributed by S. Har-Peled
# under MIT License

# Not currently used by anything...

function  print_int_w_commas( n::Int64 )
    if  ( abs( n ) < 1000 )
        print( n )
        return
    end
    print_int_w_commas( floor( Int64, n/1000 ) );
    @printf( ",%03d",  n % 1000 )
end


### Functions try to deal with floating point issues.


function  fp_ratio( a::Float64, b::Float64 )::Float64
    if   a == b
        return  1.0;
    end
    return  abs( a - b ) / (abs( a)  + abs(b) );
end

function  fp_equal( a, b )::Bool
    return  a == b;
end

function  fp_equal( a::Float64, b::Float64 )::Bool
    if   a == b
        return  true;
    end
    return  abs( a - b ) <= (0.000001* (abs( a)  + abs(b) ))
end

function  fp_equal( p::Point{D,T}, q::Point{D,T} )::Bool where {D,T}
    for  i in 1:D
        if  ! fp_equal( p[ i ], q[ i ] )
            return  false;
        end
    end
    return  true;
end

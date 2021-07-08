

#=
Programita script que obtiene las funciones b y b tilde dado un radio efectivo
=#


if length(ARGS)==0
    println("Voy a usar el parámetro radio=5.0")
else
    radio=parse(Float64, ARGS[1])
    println("Voy a usar el radio=$radio")    
end

using Cubature

sigma=1.0
h=8+1/3
#limite de integracion para el modelo gaussiano
limint=22*radio
toldura=0.0001 #mas tolerancia y la integral dura truena feo en (1,3)

function distancia2d(x1::Number,y1::Number,x2::Number,y2::Number)
    result=sqrt((x1-x2)^2+(y1-y2)^2)
    return result
end

function distancia2d(x1::Array,x2::Array)
    result=distancia2d(x1[1],x1[2],x2[1],x2[2])
    return result
end


function duro2D(x1::Number,y1::Number)
    if distancia2d(x1,y1,0,0)<=radio
        return 1.0
    else
        return 0
    end
end
    
function suave2D(x1::Number,y1::Number)
    result=exp(-distancia2d(x1,y1,0,0)/(2*radio))
    return result
end
   

function bKernelDuro2D(x1::Number,y1::Number, x2::Number,y2::Number)
    #=
    Acuerdate: Estos estan centrados en el origen, despues los vamos trasla<dando 
    sea necesario 
    =#
    dist=distancia2d(x1,y1,x2,y2)
    result=asinh(h/dist)*duro2D(x2,y2)
    return result
end

function bKernelSuave2D(x1::Number,y1::Number, x2::Number,y2::Number)
    #=
    Acuerdate: Estos estan centrados en el origen, despues los vamos trasla<dando 
    sea necesario 
    =#
    dist=distancia2d(x1,y1,x2,y2)
    result=asinh(h/dist)*suave2D(x2,y2)
    return result
end


println("funciones definidas")

xxs=-63:63
yys=-63:63

btildesuavecero=zeros(127,127)
btildeduracero=zeros(127,127) 


for j in 1:127, k in 1:127
    
    x=xxs[j]
    y=yys[k]

    btildesuavecero[k,j]=suave2D(x,y)
    btildeduracero[k,j]=duro2D(x,y)
    
end

println("funciones b tilde calculadas y guardadas en archivos")
writedlm("BtildeceroSuave-r-$radio.dat",btildesuavecero)
writedlm("BtildeceroDura-r-$radio.dat",btildeduracero)


hcubature(x ->
    bKernelDuro2D(0.001,0.001, x[1],x[2]),
[-radio,-radio],[radio,radio],reltol=toldura)


bcerodura=zeros(127,127)
bcerosuave=zeros(127,127);

for j in 64:127, k in 64:j

    (zz1,pupu)=hcubature(x ->
    bKernelDuro2D(xxs[j],yys[k], x[1],x[2]),
[-radio,-radio],[radio,radio],reltol=toldura)
    bcerodura[k,j]=zz1

    (zz2,pupu)=hcubature(x ->
    bKernelSuave2D(xxs[j],yys[k], x[1],x[2]),
[-limint,-limint],[limint,limint],reltol=0.00001)
    bcerosuave[k,j]=zz2

end


j=64
k=64
deltita=0.0000001

(bcerodura[k,j], pupu)=hcubature(x ->
bKernelDuro2D(xxs[j]+deltita,yys[k]+deltita, x[1],x[2]),
[-radio,-radio],[radio,radio],reltol=toldura)

(bcerosuave[k,j],pupu)=hcubature(x ->
    bKernelSuave2D(xxs[j]+deltita,yys[k]+deltita, x[1],x[2]),
                                 [-limint,-limint],[radio,radio],reltol=0.00001)


for j in 64:127, k in 64:j
    bcerodura[j,k]=bcerodura[k,j]
    bcerosuave[j,k]=bcerosuave[k,j]
end

for j in 64:127, k in 64:127
    bcerodura[127-k+1,127-j+1]=bcerodura[k,j]
    bcerodura[127-k+1,j]=bcerodura[k,j]
    bcerodura[k,127-j+1]=bcerodura[k,j]
    
    bcerosuave[127-k+1,127-j+1]=bcerosuave[k,j]
    bcerosuave[127-k+1,j]=bcerosuave[k,j]
    bcerosuave[k,127-j+1]=bcerosuave[k,j]
end

println("funciones b calculadas a partir de la integral eléctrica y guardadas")

writedlm("BceroDura-r-$radio.dat", bcerodura)
writedlm("BceroSuave-r-$radio.dat", bcerosuave)

using JLD

numarg=length(ARGS)
if numarg<=1
    error("Dame un nombre de archivo de funciones b y btilde para trabajar")
elseif numarg>1
    bnombre=ARGS[1] 
    btildenombre=ARGS[2]
end

if numarg==3
    archivo=ARGS[3]
    saturados=load(archivo)["CanalesSaturados"]
    push!(saturados, [1,1])
else
    saturados=Set{Array{Int,1}}()
    push!(saturados, [1,1])
end

    

println("Numero de procesos en paralelo ",nprocs())

println("Cargando las tablas de b_0 y  btilde_0...") 
b=readdlm(bnombre)
btilde=readdlm(btildenombre)

println("Iniciando las coordenadas de interes...")       
todaslasX=Array[]

#LasXNetas=Array[]
#CasiXNetas=Array[]


#alto=64
#ancho=64


alto=64
ancho=64
for j=1:alto,k=1:ancho
    push!(todaslasX,[j,k])
end


xpurgadas=filter(q->!(q in saturados), todaslasX)

numsaturados=length(saturados)
numxpurgadas=length(xpurgadas)
println("Tengo ",numsaturados, " canales saturados")
println("Tengo ",numxpurgadas, " canales para trabajar")

#=
for j=1:24, k=1:24
    push!(LasXNetas,[j+5,k+23])
end
=#

LasXchiquitas=xpurgadas[1:24]

println("Definiendo las funciones necesarias...")

@everywhere function sumaalgunasB(LasB::Array, lasX::Array,rango)
    tantas=length(rango)
    tamanho=length(lasX)
    result=zeros(tamanho,tamanho)
    xconstante=[64,64]
    for j=1:tamanho
     xj=lasX[j]
         for k=1:j
            xk=lasX[k]
            aux=0
            for l in rango
                xl=lasX[l]
                aux += LasB[(xk-xl+xconstante)...]*LasB[(xj-xl+xconstante)...]
          end
            result[j,k]=aux
        end
    end
    return result
    
end

@everywhere function sumaalgunasByBtilde(LasB::Array, LasBtilde::Array, lasX::Array,rango)
    tantas=length(rango)
    tamanho=length(lasX)
    result=zeros(tamanho,tamanho)
    xconstante=[64,64]
    for j=1:tamanho
     xj=lasX[j]
        for k=1:tamanho
            xk=lasX[k]
            aux=0
            for l in rango
                xl=lasX[l]
                aux += LasBtilde[(xk-xl+xconstante)...]*LasB[(xj-xl+xconstante)...]
          end
            result[j,k]=aux
  #  println("hola loco ", KdurasParalel[j,k])
        end
    end
    return result
    
end

@everywhere function divideArray(lasX::Array)
    #el proceso myid()==nprocs() va a ser el director de orquesta
    procesos=nprocs()
    largo=length(lasX)
    idx=myid()
    if idx==1 #nprocs()
        return 1:0
    else
        splits=[round(Int, s) for s in linspace(0,largo,procesos)]
        return splits[idx-1]+1:splits[idx]
    end
end

@everywhere function obtenKParcial(LasB,lasX)
    ranguito=divideArray(lasX)
    @show (ranguito)
    result= sumaalgunasB(LasB, lasX,ranguito)
    return result
end

@everywhere function obtenKtildeParcial(LasB,LasBTilde, lasX)
    ranguito=divideArray(lasX)
    @show (ranguito)
    result= sumaalgunasByBtilde(LasB, LasBTilde, lasX,ranguito)
    return result
end




function sacalaKenTrozos(LasB::Array, lasX::Array)
    largo=length(lasX)
    cachos=nprocs()
    aux=SharedArray{Float64}(largo,largo,cachos)
     @sync begin
        for p in procs()
            @async  aux[:,:,p]=remotecall_fetch( obtenKParcial, p,  LasB, lasX)
        end
        end #cierra el syncbegin
    result=zeros(largo,largo)
    result=sum(aux,3)
    return result
    #return result
end

function KtildeenTrozos(LasB::Array, LasBtilde::Array, lasX::Array)
    largo=length(lasX)
    cachos=nprocs()
    aux=SharedArray{Float64}(largo,largo,cachos)
     @sync begin
        for p in procs()
            @async  aux[:,:,p]=remotecall_fetch(obtenKtildeParcial, p, LasB, LasBtilde, lasX)
        end
        end #cierra el syncbegin
    result=zeros(largo,largo)
    result=sum(aux,3)
    return result
    #return result
end


println("Precompilando justo a tiempo las funciones...")

laKParcial=obtenKParcial(b,LasXchiquitas)
laKParcial=obtenKtildeParcial(b, btilde,LasXchiquitas)

K=sacalaKenTrozos(b,LasXchiquitas);
KTilde=KtildeenTrozos(b, btilde, LasXchiquitas);


println("Calculando K y Ktilde...")
@time K=sacalaKenTrozos(b,xpurgadas);
@time KTilde=KtildeenTrozos(b,btilde, xpurgadas);




println("Simetrizando K y arreglando K tilde...")


longus=size(K,1)

K=reshape(K,longus,longus)
KTilde=reshape(KTilde,longus,longus)

K=K+transpose(K)

#estabilización de matrix
#=
for j=1:longus
    K[j,j]=K[j,j]/2+0.001
end
=#


KTT=transpose(KTilde)*inv(K)

if numarg<3
    println("Escribiendo K.dat y KTilde.dat en el disco...")
    writedlm("Krevisarestable.dat",K)
    writedlm("KTilderevisarestable.dat",KTilde)

else
    nota=string("Usamos las siguientes B y BT para obtener las Ks")
    nota=string(nota,": ",ARGS[1], " y ", ARGS[2])
    println("voy a poner esta nota")
    println(nota)

    println("Escribiendo K.dat y KTilde.dat en el disco...")
    writedlm("K.dat",K)
    writedlm("KTilde.dat",KTilde)

    println("voy a modififar tu archivo ", archivo)
    println("va a tener una entrada K y KTilde")
    paguardar=load(archivo)
    paguardar["KTilde"]=KTilde
    paguardar["K"]=K
    paguardar["Nota"]=nota
    save(archivo,paguardar)
end


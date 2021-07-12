using JLD

numarg=length(ARGS)
if numarg<1
    error("Dame un nombre de archivo JLD que tenga LFP y KTT_KInv guardados")
else
    nombre=ARGS[1] 
end

println("cargando lfp")
lfp=load(nombre)["LFPSaturados"]
K=load(nombre)["K"]
KTilde=load(nombre)["KTilde"]
saturados=load(nombre)["CanalesSaturados"]
tmax=size(lfp,3)



todaslasX=Array[]


for j=1:64,k=1:64
    push!(todaslasX,[j,k])
end

xpurgadas=filter(q->!(q in saturados), todaslasX)
nbuenas=length(xpurgadas)

lambda=0.05

for j=1:nbuenas
    K[j,j]=K[j,j]/2.*(1+lambda)
end
println("K Estabilizada con lambda=",lambda)


KTT=transpose(KTilde)*inv(K)

println("KTT calculada")

CSDtentativa=zeros(4096)
CSD=zeros(lfp)

lfpv=zeros(nbuenas,tmax)

println("Acomodando los LFP correctos")

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    lfpv[j,:]=lfp[renglon,columna,:]
end


CSDTentativa=zeros(lfpv)

println("Empezando calculo")
for t=1:tmax
    CSDTentativa[:,t]=KTT*lfpv[:,t] 
end

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    CSD[renglon,columna,:]=CSDTentativa[j,:]
end


println("terminando calculo")

paguardar=load(nombre)
paguardar["kCSDCorrecta"]=CSD
save(nombre,paguardar)

println("Tu archivo jld ha sido modificado, checa una entrada kCSD")





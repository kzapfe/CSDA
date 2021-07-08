#= Version 02 porque dejaste la otra en el Labo
y se te olvido agregarla al github al parecer =#

using JLD

numarg=length(ARGS)
if numarg<2
    error("Dame un nombre de archivo JLD que tenga un KTT_KInv guardados y otro
          para calcularle su CSD")
else
    nom1=ARGS[1]
    nom2=ARGS[2]
end

KTT=load(nom1)["KTT_KInv"]
saturados=load(nom1)["CanalesSaturados"]
lfp=load(nom2)["LFPTotal"]


tmax=size(lfp,3)


todaslasX=Array[]

for j=1:64,k=1:64
    push!(todaslasX,[j,k])
end

xpurgadas=filter(q->!(q in saturados), todaslasX)
nbuenas=length(xpurgadas)

lfpv=zeros(nbuenas,tmax)

println("Acomodando los LFP correctos")

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    lfpv[j,:]=lfp[renglon,columna,:]
end


CSD=zeros(lfp)
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

paguardar=load(nom2)
paguardar["kCSDPrestado"]=CSD
save(nom2,paguardar)

println("Tu archivo jld ha sido modificado, checa una entrada kCSD")



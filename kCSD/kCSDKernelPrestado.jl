using JLD

numarg=length(ARGS)
if numarg<2
    error("Dame un nombre de archivo JLD que tenga LFP y KTT_KInv guardados y otro para trabajar")
else
    nom1=ARGS[1]
    nom2=ARGS[2]
end


datosconkernel=load(nom1)
otrosdatos=load(nom2)

#Carganto todo
lfp1=datosconkernel["LFPSaturados"];
lfp2=otrosdatos["LFPTotal"];
KTT=datosconkernel["KTT_KInv"] # Si es el mismo experimento deberian de ser los mismos.
saturados=datosconkernel["CanalesSaturados"] # Si es el mismo experimento deberian de ser los mismos.
tmax=size(lfp2,3)

#hagamos el array de electrodos validos
todaslasX=Array[]

for j=1:64,k=1:64
    push!(todaslasX,[j,k])
end

xpurgadas=filter(q->!(q in saturados), todaslasX)
nbuenas=length(xpurgadas)

#preparamos el array para obtener el CSD
CSDtentativa=zeros(4096)
CSD=zeros(lfp2)
lfp2v=zeros(nbuenas,tmax);
CSDTentativa=zeros(lfp2v);

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    lfp2v[j,:]=lfp2[renglon,columna,:]
end


for t=1:tmax
    CSDTentativa[:,t]=KTT*lfp2v[:,t] 
end

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    CSD[renglon,columna,:]=CSDTentativa[j,:]
end

#y listo

paguardar=load(nom2)
paguardar["kCSDconKernelPrestado"]=CSD
save(nom2,paguardar)

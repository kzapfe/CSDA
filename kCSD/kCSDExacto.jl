#=
Obtains the CSD from a previous calculated kernel. Needs both LFP,
and the product of both  kernels described in the paper from Potworowski.
=# 

using JLD

numarg=length(ARGS)
if numarg<1
    error("Give me a JLD that has LFP and KTT_KInv inside")
else
    nombre=ARGS[1] 
end

println("loading lfp")
lfp=load(nombre)["LFPSaturados"]
KTT=load(nombre)["KTT_KInv"]
saturados=load(nombre)["CanalesSaturados"]
tmax=size(lfp,3)

todaslasX=Array[]

for j=1:64,k=1:64
    push!(todaslasX,[j,k])
end

xpurgadas=filter(q->!(q in saturados), todaslasX)
nbuenas=length(xpurgadas)

CSDtentativa=zeros(4096)
CSD=zeros(lfp)

lfpv=zeros(nbuenas,tmax)

println("Setting up the LFP of usuable electrodes")

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    lfpv[j,:]=lfp[renglon,columna,:]
end


CSDTentativa=zeros(lfpv)

println("Convoluting...")
for t=1:tmax
    CSDTentativa[:,t]=KTT*lfpv[:,t] 
end

for j=1:nbuenas
    renglon=xpurgadas[j][1]
    columna=xpurgadas[j][2]
    CSD[renglon,columna,:]=CSDTentativa[j,:]
end


println("Finishing  operations")

paguardar=load(nombre)
paguardar["kCSDCorrecta"]=CSD
save(nombre,paguardar)

println("Your file has been modified. Look for kCSD entry.")





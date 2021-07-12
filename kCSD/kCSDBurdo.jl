using JLD

numarg=length(ARGS)
if numarg<1
    error("Dame un nombre de archivo JLD que tenga LFP guardados")
else
    nombre=ARGS[1] 
end

println("cargando lfp")
lfp=load(nombre)["LFPSaturados"]

tmax=size(lfp,3)

if numarg==1
 if isfile("KTT_Suave_r_0.5_4095.dat")
     KTT=readdlm("KTT_Suave_r_0.5_4095.dat")
     println("KTT leida del disco")
 else
     K=readdlm("K_Suave_r-0.5-4095.dat")
     KTilde=readdlm("KTilde_Suave_r-0.5-4095.dat")
     Kinv=inv(K)
     KTT=transpose(KTilde)*Kinv;
     writedlm("KTT_Suave_r_0.5_4095.dat", KTT)
     println("KTT Kinv calculada y guardada")
 end
    nombrekcsd="kCSD_Suave_r_0_5"
elseif numarg==2
    
    println("voy a usar una KTT del disco duro que tu dijiste")
    palabra=ARGS[2]

    KTT=readdlm(palabra)
    
    nombrekcsd=replace(palabra, "KTT", "kCSD")
    nombrekcsd=replace(nombrekcsd, ".dat", "")
else

    println("voy a usar K y KT del disco duro que tu dijiste")
    
    palabra=ARGS[2]
    
    K=readdlm(ARGS[2])
    KTilde=readdlm(ARGS[3])
    Kinv=inv(K)
    KTT=transpose(KTilde)*Kinv;   
   
    nombrekcsd=replace(palabra, "K_", "kCSD_")
    nombrekcsd=replace(nombrekcsd, ".dat", "")

end



CSDtentativa=zeros(4096)
CSD=zeros(lfp)

println("Empezando calculo")
for t=1:tmax
    LFPvector1=reshape(transpose(lfp[:,:,t]),4096);
    lfpvecredux=LFPvector1[2:end]
    CSDtentativa[2:end]=transpose(KTT)*lfpvecredux
    CSD[:,:,t]=reshape(CSDtentativa,64,64)
end

println("terminando calculo")

paguardar=load(nombre)
paguardar[nombrekcsd]=CSD
save(nombre,paguardar)

println("Tu archivo jld ha sido modificado, checa una entrada kCSD")





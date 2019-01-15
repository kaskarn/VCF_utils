using CodecZlib, GeneticVariation, DataFrames, CSV, GLM, Distributions

vcfpath="/proj/epi/Genetic_Data_Center/calico/PAGEII_genotypes/MEGA/imputed_data/vcf/chr01.vcf.gz"
idpath="/proj/epi/CVDGeneNas/antoine/CHARGE_ECG_GWAS/phenotypes/hchs_unrelated_ids.txt"
outpath="/proj/epi/CVDGeneNas/antoine/CHARGE_ECG_GWAS/output/hwe_res.txt"

const vcfind = idrdy[:ind]
const chi2 = Chisq(1)

function hwe_var(vcfnow, vcfind)
    gparr = zeros(Float64, 3)

    gnow = VCF.genotype(vcfnow, vcfind, "GP");
    for i in eachindex(gnow)
        gparr += parse.(Float64, split(gnow[i],','))
    end
    NAA, NAa, Naa = gparr
    chi2_score = length(vcfind) * ((4*NAA*Naa - (NAa)^2) / ((2NAA + NAa)*(2Naa+NAa)))^2
    p = 1 - cdf(chi2, chi2_score)
    VCF.id(vcfnow)[1], NAA, NAa, Naa, p
end

function hwe_vcf(reader, out, vcfind)
    for vcfnow in reader
        tmp = hwe_var(vcfnow, vcfind)
        join(out, tmp, '\t')
        write(out, '\n')
        println(tmp[1])
    end
end

reader = VCF.Reader(GzipDecompressorStream(open(vcfpath)));
ids = DataFrame(readdlm(idpath, String))
rename!(ids, :x1 => :id_rdy)

vcfids = DataFrame(id_rdy = reader.header.sampleID, ind = 1:length(reader.header.sampleID))
idrdy = join(ids, vcfids, on = :id_rdy, kind = :inner)

out = open(outpath, "w")
write(out, "snp\tN_AA\tN_Aa\tN_aa\tp_hwe\n")
hwe_vcf(reader, out, vcfind)

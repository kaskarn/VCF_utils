using CodecZlib, GeneticVariation, DataFrames, CSV, GLM

#parsing command line
incmd = join(ARGS, " ")
pinputs = Dict()
[get!(pinputs, split(i)[1], size(split(i),1) > 1 ? split(i)[2] : true) for i in split(incmd, "--")[2:end]]

[println("$(k): $(pinputs[k])") for k in keys(pinputs)];

reader = haskey(pinputs, "vcfpath") ? VCF.Reader(GzipDecompressorStream(open(pinputs["vcfpath"]))) : VCF.Reader(open(STDIN))
pheno = CSV.read(pinputs["phepath"], delim = '\t', null = "NA", rows_for_type_detect = 10000) #SUBOPTIMAL UNTIL UPDATE TO CSV PARSER :(

vcfids = DataFrame(id_rdy = reader.header.sampleID, ind = 1:length(reader.header.sampleID))
rename!(pheno, Symbol(pinputs["phenid"]) => :id_rdy)
vcf_phen = join(pheno, vcfids, on = :id_rdy, kind = :inner)
sort!(vcf_phen , cols = order(:ind))
vcfind = vcf_phen[:ind]

formula = Formula(Symbol(pinputs["phen"]), parse(join(["G + ", pinputs["rhs"]])))

outres=STDOUT
(pinputs["outfile"] == "STDOUT") || (outres = open(pinputs["outfile"], "w+"))

write(outres, "VCF_ID\tCHROM\tPOS\tREF\tALT\tALT_AF\tALT_AC\tN_INFORMATIVE\tBETA\tSE\tPVALUE")

gind = trues(length(vcfind))
vcf_phen[:G] = rand(size(vcf_phen,1))
tmdf = fit(LinearModel, formula, vcf_phen).mf.df
mf = convert(Matrix{Float64}, Array(tmdf))
mf[:,1] = 1
y = vcf_phen[Symbol(pinputs["phen"])]

function process_var!(vcfnow, y, mf, gind, vcfind, outres)
	gnow = tryparse.(Float64, VCF.genotype(vcfnow, vcfind, "DS"))
	n, ac = 0, 0.0
	for (i, g) in enumerate(gnow)
		gind[i] = g.hasvalue
		if g.hasvalue
			n += 1
			ac += g.value
			mf[i,2] = g.value
		end
	end

	write(outres, '\n')
	(ac < 2.0 || ac > 2n-2) && return join(outres, [VCF.id(vcfnow)[1], VCF.chrom(vcfnow), VCF.pos(vcfnow), VCF.ref(vcfnow), VCF.alt(vcfnow)[1], ac/2/n, round(ac), n, "NA", "NA", "NA"], '\t')
	mnow = lm(mf[gind,:],y[gind])
	b = coef(mnow)[2]
	se = stderr(mnow)[2]
	p = GLM.ccdf(GLM.FDist(1,GLM.df_residual(mnow)),abs2(b/se))
	return join(outres, [VCF.id(vcfnow)[1], VCF.chrom(vcfnow), VCF.pos(vcfnow), VCF.ref(vcfnow), VCF.alt(vcfnow)[1], ac/2/n, round(ac), n, b, se, p], '\t')
end

println("GWAS Started on $(size(vcf_phen, 1)) Individuals\n\n")

i = 0
haskey(pinputs, "test") && (pinputs["test"] = parse(Int64, pinputs["test"]))
!haskey(pinputs, "test") && (get!(pinputs, "test", typemax(Int64)))
@time for vcfnow in Iterators.take(reader, pinputs["test"])
	i += 1
	process_var!(vcfnow, y, mf, gind, vcfind, outres)
	i%10000 == 0 && println("$(i) variants processed")
end

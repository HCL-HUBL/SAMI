process splicing_nosplice {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	tuple path("I.rds"), path("S.rds"), path("groups.rds"), path("sites.rds"), path("events.rds")
	path("depth.bed")
	
	output:
	tuple path("out/I.rds"), path("out/S.rds"), path("out/groups.rds"), path("out/sites.rds"), path("out/events.rds"), emit: RDS

	shell:
	template 'splicing_nosplice.R'
}

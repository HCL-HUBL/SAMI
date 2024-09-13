process splicing_nosplice {
	cpus 1
	time { 10.minute * task.attempt }
	memory { 1.GB * task.attempt }

	input:
	tuple path("I.rds"), path("S.rds"), path("groups.rds"), path("sites.rds"), path("events.rds")
	path("depth/*")
	
	output:
	tuple path("out/I.rds"), path("out/S.rds"), path("out/groups.rds"), path("out/sites.rds"), path("out/events.rds"), path("out/depth.rds"), emit: RDS

	shell:
	template 'splicing_nosplice.R'
}

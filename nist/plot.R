
library(plotly)
data = read.csv('nist.csv', as.is=T)[-1,]
data$T = as.numeric(data$T)
data$p = as.numeric(data$p)
data$V = as.numeric(data$V)
data$k = as.numeric(data$k)
data$E = as.numeric(data$E)
data$H = as.numeric(data$H)
data$S = as.numeric(data$S)
data$Cv = as.numeric(data$Cv)
data$Cp = as.numeric(data$Cp)
data$u = as.numeric(data$u)
data$JT = as.numeric(data$JT)
data$nu = as.numeric(data$nu)
data$rho = as.numeric(data$rho)

gas = data[data$phase=='vapor',]
liquid = data[data$phase=='liquid',]
supercritical = data[data$phase=='supercritical',]

#namespace
Evaluation = function(observations, property, predict){
	input = observations
	output = observations[,property]
	list(
		pe = function(parameters){
			(output-predict(input, parameters))/output
		},
		e = function(parameters){
			output-predict(input, parameters)
		},
		plot = function(parameters1, parameters2) {
			plot_ly(
				x=c(input$T, input$T, input$T), 
				y=log(c(input$p, input$p, input$p)), 
				z=c(output, predict(input, parameters1), predict(input, parameters2)), 
				color=c(
					rep('observed', times=length(output)), 
					rep('parameters1', times=length(output)), 
					rep('parameters2', times=length(output))
				)
			) %>% layout(scene = list(xaxis=list(title='T'), yaxis=list(title='p')))
		}
	)
}

#namespace
Costs = function(evaluation){
	list(
		# mean percent error
		mpe = function(parameters) {
			mean(abs(evaluation$pe(parameters)))
		},
		# mean error
		me = function(parameters) {
			mean(abs(evaluation$e(parameters)))
		},
		# max percent error
		mxpe = function(parameters) {
			max(abs(evaluation$pe(parameters)))
		},
		# max error
		mxe = function(parameters) {
			max(abs(evaluation$e(parameters)))
		}
	)
}

Iteration = function(observations, evaluation, costs, model, cost_id){
	list(
		initialize = function(){ list( par=model$initialize(), value=costs[[cost_id]](model$initialize()) ) },
		iterate = function(state){ optim( state$par, costs[[cost_id]] ) },
		summary = function(state){ 
			start = model$initialize()
			print(state$par)
			print(paste(model$code(state$par), ' \r\n// ', 
				unique(observations$compound), ', ',
				'mean error: ', costs$mpe(state$par), 
				'max error: ', costs$mxpe(state$par), 
				'range: ', min(observations$T, na.rm=TRUE), '-', max(observations$T, na.rm=TRUE), 'K, ', 
				           min(observations$p, na.rm=TRUE), '-', max(observations$p, na.rm=TRUE), 'MPa, ', 
                'stp estimate: ', sprintf('%.3f', model$predict(list(p=0.1, T=273.15), state$par)), 
                sep=''))
			evaluation$plot(start, state$par)
		}
	)
}

Fitting = function(property, model, cost_id){
	list(
		optimize = function(observations){
			evaluation = Evaluation(observations, property, model$predict)
			costs = Costs(evaluation)
			iteration = Iteration(observations, evaluation, costs, model, cost_id)
			iteration$summary(
				iteration$iterate(
					iteration$initialize()))
		},
		plot = function(observations){
			plot_ly(
				x=observations$T, 
				y=log(observations$p), 
				z=observations[,property], 
				color=observations$compound
			) %>% layout(scene = list(xaxis=list(title='T'), yaxis=list(title='p')))
		}
	)
}


gas_heat_capacity_hydrogen = list(
	initialize=function(){ c(1e-3,1,  1e-6,1.7,  8,120,100,  6) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		x7 = parameters[7]
		x8 = parameters[8]
		p = data$p
		T = data$T
		u = (T-x6)/x7
		return( x1*p^x2 + x3*T^x4 + x5*u/sqrt(1+u^2) + x8 )
	},
	code = function(parameters){
		paste('get_sigmoid_exponent_pressure_temperature_function\r\n',
			      '(si::kelvin, si::megapascal, si::joule/(si::gram * si::kelvin),\r\n', 
			      paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

gas_heat_capacity_dipper102 = list(
	initialize=function(){ c(2.7, 0.35, -1.8, 2.6,    1e-3, 0.75) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		p = data$p
		T = data$T
		return( x1*T^x2 / (1+x3/T + x4/T^2)  + x5*T^x6 )
	},
	code = function(parameters){
		paste('get_dippr_temperature_pressure_function_102\r\n',
			'(si::kelvin, si::megapascal, si::joule/(si::gram * si::kelvin),\r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

hydrogen = list(
	predict=gas_heat_capacity_dipper102$predict,
	code=gas_heat_capacity_dipper102$code,
	initialize=function(){ c(0.1,1.02,-2.26,1.56,    -0.181, 0.640) }
)
gas_heat_capacities = Fitting('Cp', gas_heat_capacity_hydrogen, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'hydrogen'           ,])




# For a given compound and pressure, filter out all but the measurements taken above a temperature at which the lowest measurement occurs.
gas.Cp.id.df = aggregate(
  seq(nrow(gas)),
  list(gas$phase, gas$compound, gas$p), 
  function(ids) {  
    subset = gas[ids,]
    ids[ subset$T > subset$T[  subset$Cp == min(subset$Cp)  ] ]   
  }
)
gas.Cp = gas[stack(gas.Cp.id.df$x)$values,] 

gas_heat_capacity_model1 = list(
	initialize=function(){ c(1e-2,1,  1e-4,0.8,  1) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*T^x4  + x5 )
	},
	code = function(parameters){
		paste('get_exponent_pressure_temperature_function\r\n',
			'(si::kelvin, si::megapascal, si::joule/(si::gram * si::kelvin),\r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

gas_heat_capacity_sigmoid = list(
	initialize=function(){ c(1e-3,1,  1,160,1e-3,  0) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*((T-x5)/x4)/sqrt(1+((T-x5)/x4)^2) + x6 )
	},
	code = function(x){
		paste('get_sigmoid_exponent_pressure_temperature_function\r\n',
			      '(si::kelvin, si::megapascal, si::joule/(si::gram * si::kelvin),\r\n', 
			      paste(sprintf('%.5f', c(x[1],x[2], 0,0, x[3],x[4],x[5], x[6])), collapse=', '), '),', sep='')
	}
)

gas_heat_capacity_sigmoid2 = list(
	initialize=function(){ c(1e-3,1,  1,160,400, 2) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*((T-x5)/x4)/sqrt(1+((T-x5)/x4)^2) + x6 )
	},
	code = function(x){
		paste('get_sigmoid_exponent_pressure_temperature_function\r\n',
		      '(si::kelvin, si::megapascal, si::joule/(si::gram * si::kelvin),\r\n', 
		      paste(sprintf('%.5f', c(x[1],x[2], 0,0, x[3],x[4],x[5], x[6])), collapse=', '), '),', sep='')
	}
)


gas_heat_capacity_hydrogen = list(
	initialize=function(){ c(1e-3,1,  0.001,1.7,  8,120,100,  6) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		x7 = parameters[7]
		x8 = parameters[8]
		p = data$p
		T = data$T
		return( x1*p^x2 + (T/x3)^x4 + x5*((T-x6)/x7)/sqrt(1+((T-x6)/x7)^2) + x8 )
	},
	code = function(parameters){
		paste('get_sigmoid_exponent_pressure_temperature_function\r\n',
			      '(si::kelvin, si::megapascal, si::joule/(si::gram * si::kelvin),\r\n', 
			      paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)


gas_heat_capacities$plot(gas.Cp[gas.Cp$compound == 'hydrogen',])

gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'water'              ,]) 
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'carbon monoxide'    ,]) 
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'ethane'             ,]) 
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'ammonia'            ,]) 
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'sulfur dioxide'     ,])
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'tetrafluoromethane' ,])
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'benzene'            ,])

gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'nitrogen'           ,])
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'oxygen'             ,])
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'carbon dioxide'     ,])

gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid2, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'methane'            ,])

gas_heat_capacities$optimize(data[data$compound == 'argon'              ,]) # try constant
gas_heat_capacities$optimize(data[data$compound == 'helium'             ,]) # try constant

gas_heat_capacities$plot(gas.Cp[gas.Cp$compound == 'helium',])








# For a given compound and pressure, filter out all but the measurements taken above a temperature at which the lowest measurement occurs.
gas.k.id.df = aggregate(
  seq(nrow(gas)),
  list(gas$phase, gas$compound, gas$p), 
  function(ids) {  
    subset = gas[ids,]
    ids[ subset$T > subset$T[  subset$k == min(subset$k)  ] ]   
  }
)
gas.k = gas[stack(gas.k.id.df$x)$values,]

gas_conductivity_model = list(
	initialize=function(){ c(1e-3,1.2,  2.5e-4,0.8,  0) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*T^x4  + x5 )
	},
	code = function(parameters){
		paste('get_exponent_pressure_temperature_function\r\n',
			'(si::kelvin, si::megapascal, si::watt/(si::meter * si::kelvin),\r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

gas_conductivity_sigmoid1 = list(
	initialize=function(){ c(1e-3,1.1,  0.07,600,1200, 0.07) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*((T-x5)/x4)/sqrt(1+((T-x5)/x4)^2) + x6 )
	},
	code = function(x){
		paste('get_sigmoid_exponent_pressure_temperature_function\r\n',
		      '(si::kelvin, si::megapascal, si::watt/(si::meter * si::kelvin),\r\n', 
		      paste(sprintf('%.5f', c(x[1],x[2], 0,0, x[3],x[4],x[5], x[6])), collapse=', '), '),', sep='')
	}
)

gas_conductivity_sigmoid2 = list(
	initialize=function(){ c(1e-3, 1.1, 0.07, 500, 500, 0.07) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		x6 = parameters[6]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*((T-x5)/x4)/sqrt(1+((T-x5)/x4)^2) + x6 )
	},
	code = function(x){
		paste('get_sigmoid_exponent_pressure_temperature_function\r\n',
		      '(si::kelvin, si::megapascal, si::watt/(si::meter * si::kelvin),\r\n', 
		      paste(sprintf('%.5f', c(x[1],x[2], 0,0, x[3],x[4],x[5], x[6])), collapse=', '), '),', sep='')
	}
)

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'carbon dioxide'    ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'nitrogen'          ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'oxygen'            ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'argon'             ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'hydrogen'          ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'helium'            ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'carbon monoxide'   ,])
gas.conductivities$optimize(gas.k[gas.k$compound == 'tetrafluoromethane',])

gas.conductivities = Fitting('k', gas_conductivity_sigmoid1, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'water'             ,]) # needs custom fit
gas.conductivities$optimize(gas.k[gas.k$compound == 'methane'           ,]) # needs custom fit

gas.conductivities = Fitting('k', gas_conductivity_sigmoid2, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'ammonia'           ,]) # needs custom fit
gas.conductivities$optimize(gas.k[gas.k$compound == 'ethane'            ,]) # needs custom fit

# gas.conductivities$optimize(gas.k[gas.k$phase=='vapor' & gas.k$compound == 'benzene'            & (gas.k$p<3 | gas.k$T>145),])
# gas.conductivities$optimize(gas.k[gas.k$phase=='vapor' & gas.k$compound == 'sulfur dioxide'     & (gas.k$p<3 | gas.k$T>145),])









gas.nu.id.df = aggregate(
  seq(nrow(gas)),
  list(gas$phase, gas$compound, gas$p), 
  function(ids) {  
    subset = gas[ids,]
    ids[ subset$T > subset$T[  subset$nu == min(subset$nu)  ] ]   
  }
)
gas.nu = gas[stack(gas.nu.id.df$x)$values,]


gas_viscosity_model = list(
	initialize=function(){ c(0.1,1.5,  0.3,0.8,  0) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*T^x4  + x5 )
	},
	code = function(parameters){
		paste('get_exponent_pressure_temperature_function\r\n',
			'(si::kelvin, si::megapascal, si::micropascal*si::second, \r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'water'              ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'carbon dioxide'     ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'nitrogen'           ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'oxygen'             ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'argon'              ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'methane'            ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'helium'             ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'hydrogen'           ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'carbon monoxide'    ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'ethane'             ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'ammonia'            ,])
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'tetrafluoromethane' ,])

# gas_viscosities$optimize(gas.nu[gas.nu$phase=='vapor' & gas.nu$compound == 'benzene'            & (gas.nu$p<3 | gas.nu$T>145),])
# gas_viscosities$optimize(gas.nu[gas.nu$phase=='vapor' & gas.nu$compound == 'sulfur dioxide'     & (gas.nu$p<3 | gas.nu$T>145),])






#reorganize to make copy pasting results into C++ easier

gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'water'              ,])

gas.conductivities = Fitting('k', gas_conductivity_sigmoid1, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'water'             ,]) # needs custom fit

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'water'              ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'nitrogen'           ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'nitrogen'          ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'nitrogen'           ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'oxygen'             ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'oxygen'            ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'oxygen'             ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'carbon dioxide'     ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'carbon dioxide'    ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'carbon dioxide'     ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid2, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'methane'            ,])

gas.conductivities = Fitting('k', gas_conductivity_sigmoid1, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'methane'           ,]) # needs custom fit

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'methane'            ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid2, 'mpe')
gas_heat_capacities$optimize(data[data$compound == 'argon'              ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'argon'             ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'argon'              ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_sigmoid2, 'mpe')
gas_heat_capacities$optimize(data[data$compound == 'helium'             ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'helium'            ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'helium'             ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_hydrogen, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'hydrogen'           ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'hydrogen'          ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'hydrogen'           ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'ammonia'            ,])

gas.conductivities = Fitting('k', gas_conductivity_sigmoid2, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'ammonia'           ,]) # needs custom fit

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'ammonia'            ,])



gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'sulfur dioxide'     ,])




gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'carbon monoxide'    ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'carbon monoxide'   ,])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'carbon monoxide'    ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'ethane'             ,])

gas.conductivities = Fitting('k', gas_conductivity_sigmoid2, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'ethane'            ,]) # needs custom fit

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'ethane'             ,])


gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'tetrafluoromethane' ,])

gas.conductivities = Fitting('k', gas_conductivity_model, 'mpe')
gas.conductivities$optimize(gas.k[gas.k$compound == 'tetrafluoromethane',])

gas_viscosities = Fitting('nu', gas_viscosity_model, 'mpe')
gas_viscosities$optimize(gas.nu[gas.nu$compound == 'tetrafluoromethane' ,])






gas_heat_capacities = Fitting('Cp', gas_heat_capacity_model1, 'mpe')
gas_heat_capacities$optimize(gas.Cp[gas.Cp$compound == 'benzene'            ,])






gas_densities = Fitting('rho', gas_viscosity_model, 'mpe')
gas_densities$plot(gas)














liquid_heat_capacity_model1 = list(
	initialize=function(){ c(1e-3,1.0,  1e-3,1.0,  1.0) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		p = data$p
		T = data$T
		return( x1*p^x2  + x3*T^x4  + x5 )
	},
	code = function(parameters){
		paste('get_exponent_pressure_temperature_function\r\n',
			'(si::kelvin, si::megapascal, si::watt/(si::meter * si::kelvin), \r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

liquid_heat_capacity_model2 = list(
	initialize=function(){ c(1e-1,1000.0,  0.5,300.0,  0.5) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		x5 = parameters[5]
		p = data$p
		T = data$T
		return( x1*exp(p/x2) + x3*exp(T/x4)  + x5 )
	},
	code = function(parameters){
		paste('get_exponential_pressure_temperature_function\r\n',
			'(si::kelvin, si::megapascal, si::watt/(si::meter * si::kelvin), \r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

liquid_heat_capacity_model_perry = list(
	initialize=function(){ c(1.0, 1e-5,  1e-5, 3e-6) },
	predict=function(data, parameters){
		x1 = parameters[1]
		x2 = parameters[2]
		x3 = parameters[3]
		x4 = parameters[4]
		p = data$p
		T = data$T
		return( x1 + x2*T + x3/(T*T) + x4*T*T )
	},
	code = function(parameters){
		paste('get_exponential_pressure_temperature_function\r\n',
			'(si::kelvin, si::megapascal, si::watt/(si::meter * si::kelvin), \r\n', 
			paste(sprintf('%.5f', parameters), collapse=', '), '),', sep='')
	}
)

liquid_heat_capacities = Fitting('Cp', liquid_heat_capacity_model_perry, 'mpe')
liquid_heat_capacities$optimize(liquid[liquid$compound == 'water',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'nitrogen',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'oxygen',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'carbon dioxide',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'methane',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'argon',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'helium',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'hydrogen',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'ammonia',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'sulfur dioxide',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'carbon monoxide',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'ethane',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'tetrafluoromethane',])
liquid_heat_capacities$optimize(liquid[liquid$compound == 'benzene',])



liquid_heat_capacities = Fitting('Cp', liquid_heat_capacity_model_perry, 'mpe')

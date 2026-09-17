o2physics_add_dpl_workflow(jpsi-fragmentation-function
                    SOURCES jpsiFragmentationFunction.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(pid-efficiency
                    SOURCES PIDefficiencyConverter.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::PWGDQCore O2Physics::MLCore
                    COMPONENT_NAME Analysis)
                    
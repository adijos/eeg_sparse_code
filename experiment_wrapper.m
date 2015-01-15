function experiment_wrapper()
    datatype = 1; % 0 is preictal, 1 is interictal
    pars = default_pars(struct);
    sparseeg(datatype,pars)
end
function data = SImportMesh(data, filename)

Domain = import_bulk_model(filename);
Domain.Materials = [1 1 1.2 340];
data.Model.Domain = Domain;
data.Model = SInitModel(data.ModelType, data.Model);

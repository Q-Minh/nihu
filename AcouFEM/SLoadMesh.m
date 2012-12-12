function data = SLoadMesh(data, filename)

load(filename, 'Domain');
data.Model.Domain = Domain;
data.Model = SInitModel(data.ModelType, data.Model);

import docker
import synapseclient

parentid = "syn16925027"

syn = synapseclient.Synapse()
syn.login()

results = syn.tableQuery("SELECT id FROM syn22280982 where status = 'ACCEPTED' AND createdBy <> 3360851")
for row in results:
    submissionid = row[0]
    print(submissionid)
    sub = syn.getSubmission(submissionid)
    client = docker.from_env()
    docker_image = sub.dockerRepositoryName + "@" + sub.dockerDigest
    image = client.images.pull(docker_image)
    new_image_name = "docker.synapse.org/{}/{}".format(parentid, submissionid)
    image.tag(new_image_name)
    client.images.push(new_image_name)
    repos = list(syn.getChildren(parentid, includeTypes=['dockerrepo']))
    for repo in repos:
        name = syn.get(repo['id']).repositoryName
        if name == new_image_name:
            entity_id = repo['id']
            break







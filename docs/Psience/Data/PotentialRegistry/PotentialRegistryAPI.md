## <a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI">PotentialRegistryAPI</a> 

<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry.py#L14)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry.py#L14?message=Update%20Docs)]
</div>









<div class="collapsible-section">
 <div class="collapsible-section collapsible-section-header" markdown="1">
## <a class="collapse-link" data-toggle="collapse" href="#methods" markdown="1"> Methods and Properties</a> <a class="float-right" data-toggle="collapse" href="#methods"><i class="fa fa-chevron-down"></i></a>
 </div>
 <div class="collapsible-section collapsible-section-body collapse show" id="methods" markdown="1">
 ```python
request_base: str
potential_registry_org: str
blacklist_repos: list
```
<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.__init__" class="docs-object-method">&nbsp;</a> 
```python
__init__(self, token=None, request_delay_time=None, release_manager=None, **opts): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry.py#L17)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry.py#L17?message=Update%20Docs)]
</div>


<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.list_repos" class="docs-object-method">&nbsp;</a> 
```python
list_repos(self, owner=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L28)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L28?message=Update%20Docs)]
</div>


<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.list_potentials" class="docs-object-method">&nbsp;</a> 
```python
list_potentials(self): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L32)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L32?message=Update%20Docs)]
</div>


<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.list_releases" class="docs-object-method">&nbsp;</a> 
```python
list_releases(self, repo_or_owner, repo=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L35)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L35?message=Update%20Docs)]
</div>


<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.latest_release" class="docs-object-method">&nbsp;</a> 
```python
latest_release(self, repo_or_owner, repo=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L40)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L40?message=Update%20Docs)]
</div>


<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.get_release_list" class="docs-object-method">&nbsp;</a> 
```python
get_release_list(self, repo_or_owner, repo=None, update=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L46)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L46?message=Update%20Docs)]
</div>


<a id="Psience.Data.PotentialRegistry.PotentialRegistryAPI.get_potential" class="docs-object-method">&nbsp;</a> 
```python
get_potential(self, name, update=None, version=None): 
```
<div class="docs-source-link" markdown="1">
[[source](https://github.com/McCoyGroup/Psience/blob/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L76)/
[edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry/PotentialRegistryAPI.py#L76?message=Update%20Docs)]
</div>
 </div>
</div>












---


<div markdown="1" class="text-secondary">
<div class="container">
  <div class="row">
   <div class="col" markdown="1">
**Feedback**   
</div>
   <div class="col" markdown="1">
**Examples**   
</div>
   <div class="col" markdown="1">
**Templates**   
</div>
   <div class="col" markdown="1">
**Documentation**   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
  <div class="row">
   <div class="col" markdown="1">
[Bug](https://github.com/McCoyGroup/Psience/issues/new?title=Documentation%20Improvement%20Needed)/[Request](https://github.com/McCoyGroup/Psience/issues/new?title=Example%20Request)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/examples/Psience/Data/PotentialRegistry/PotentialRegistryAPI.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/examples/Psience/Data/PotentialRegistry/PotentialRegistryAPI.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/gh-pages/ci/docs/Psience/Data/PotentialRegistry/PotentialRegistryAPI.md)/[New](https://github.com/McCoyGroup/Psience/new/gh-pages/?filename=ci/docs/templates/Psience/Data/PotentialRegistry/PotentialRegistryAPI.md)   
</div>
   <div class="col" markdown="1">
[Edit](https://github.com/McCoyGroup/Psience/edit/master/Psience/Data/PotentialRegistry.py#L14?message=Update%20Docs)   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
   <div class="col" markdown="1">
   
</div>
</div>
</div>
</div>
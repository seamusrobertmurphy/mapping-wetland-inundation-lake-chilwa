<!DOCTYPE qgis>
<!-- Categorised style for ground truth points, keyed on the integer field
     "code". Load points as a delimited-text or GeoJSON layer, then
     Properties, Symbology, Load Style, and pick this file. -->
<qgis version="3.34" styleCategories="Symbology">
  <renderer-v2 type="categorizedSymbol" attr="code">
    <categories>
      <category render="true" symbol="0" value="1" label="open water"/>
      <category render="true" symbol="1" value="2" label="flooded vegetation"/>
      <category render="true" symbol="2" value="3" label="dry vegetation"/>
      <category render="true" symbol="3" value="4" label="bare soil"/>
      <category render="true" symbol="4" value="8" label="rice paddy"/>
      <category render="true" symbol="5" value="9" label="settlement"/>
      <category render="true" symbol="6" value="0" label="unassigned"/>
    </categories>
    <symbols>
      <symbol name="0" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="31,110,140,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
      <symbol name="1" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="95,136,54,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
      <symbol name="2" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="176,132,71,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
      <symbol name="3" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="140,109,82,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
      <symbol name="4" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="126,87,160,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
      <symbol name="5" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="138,143,148,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
      <symbol name="6" type="marker" alpha="1">
        <layer class="SimpleMarker">
          <Option type="Map">
            <Option name="color" type="QString" value="201,207,204,255"/>
            <Option name="outline_color" type="QString" value="255,255,255,255"/>
            <Option name="outline_width" type="QString" value="0.4"/>
            <Option name="size" type="QString" value="2.6"/>
            <Option name="name" type="QString" value="circle"/>
          </Option>
        </layer>
      </symbol>
    </symbols>
  </renderer-v2>
</qgis>

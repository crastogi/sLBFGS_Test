//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vJAXB 2.1.10 in JDK 6 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2015.03.06 at 02:30:00 PM EST 
//


package config;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for anonymous complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType>
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element ref="{}ExperimentReference"/>
 *       &lt;/sequence>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "", propOrder = {
    "experimentReference"
})
@XmlRootElement(name = "CrossValidation")
public class CrossValidation {

    @XmlElement(name = "ExperimentReference", required = true)
    protected ExperimentReference experimentReference;

    /**
     * Gets the value of the experimentReference property.
     * 
     * @return
     *     possible object is
     *     {@link ExperimentReference }
     *     
     */
    public ExperimentReference getExperimentReference() {
        return experimentReference;
    }

    /**
     * Sets the value of the experimentReference property.
     * 
     * @param value
     *     allowed object is
     *     {@link ExperimentReference }
     *     
     */
    public void setExperimentReference(ExperimentReference value) {
        this.experimentReference = value;
    }

}
